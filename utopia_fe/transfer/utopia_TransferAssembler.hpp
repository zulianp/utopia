#ifndef UTOPIA_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_TRANSFER_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_libmesh.hpp"
#include "utopia.hpp"
#include "utopia_Path.hpp"

#include <memory>

namespace utopia {

	class TransferOptions {
	public:
		TransferOptions()
		: from_var_num(0), to_var_num(0), n_var(1), tags({})
		{}

		int from_var_num;
		int to_var_num;
		int n_var;
		std::vector< std::pair<int, int> > tags;
	};

	class TransferAssembler final {
	public:
		using FunctionSpace = utopia::LibMeshFunctionSpace;
		using SparseMatrix  = utopia::DSMatrixd;
		using MeshBase      = libMesh::MeshBase;
		using DofMap        = libMesh::DofMap;

		class Algorithm {
		public:
			virtual ~Algorithm() {}
			// virtual bool assemble(SparseMatrix &B) = 0;
			virtual bool assemble(std::vector<std::shared_ptr<SparseMatrix> > &B) = 0;
		};

		TransferAssembler(
			const std::shared_ptr<LocalAssembler> &assembler,
			const std::shared_ptr<Local2Global>   &local2global);

		bool assemble(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			SparseMatrix &B,
			const TransferOptions &opts = TransferOptions()
		);

		bool assemble(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			std::vector<std::shared_ptr<SparseMatrix> > &B,
			const TransferOptions &opts = TransferOptions()
		);

		void set_assembler(const std::shared_ptr<LocalAssembler> &assembler)
		{
			assembler_ = assembler;
		}

		void set_local_2_global(const std::shared_ptr<Local2Global> &local2global)
		{
			local2global_ = local2global;
		}

		~TransferAssembler();

	private:
		std::shared_ptr<LocalAssembler> assembler_;
		std::shared_ptr<Local2Global> local2global_;
		std::shared_ptr<Algorithm> algorithm_;
	};

	class TransferOperator {
	public:
		virtual ~TransferOperator() {}
		virtual void apply(const DVectord &from, DVectord &to) const = 0;
		virtual void apply_transpose(const DVectord &from, DVectord &to) const = 0;
		virtual void describe(std::ostream &) const {}
		virtual bool write(const Path &) const { return false; }
	};


	/**
	 * @brief constructed as (D^-1 * B) * ( . )
	 */
	class L2TransferOperator : public TransferOperator {
	public:
		inline void apply(const DVectord &from, DVectord &to) const override
		{
			DVectord B_from = *B * from;

			if(empty(to)) {
				to = B_from;
			}

			linear_solver->apply(B_from, to);
		}

		void fix_mass_matrix_operator()
		{
			DVectord d;

			Size s = local_size(*D);
			d = local_values(s.get(0), 1.);

			{
				Write<DVectord> w_d(d);

				each_read(*D, [&d](const SizeType i, const SizeType, const double) {
					d.set(i, 0.);
				});
			}

			(*D) += DSMatrixd(diag(d));
		}

		///@brief assumes that D is symmetric
		void apply_transpose(const DVectord &from, DVectord &to) const override
		{
			DVectord D_inv_from = local_zeros(local_size(*D).get(0));
			linear_solver->apply(from, D_inv_from);
			to = transpose(*B) * D_inv_from;
		}

		inline L2TransferOperator(
			const std::shared_ptr<DSMatrixd> &B,
			const std::shared_ptr<DSMatrixd> &D,
			const std::shared_ptr<LinearSolver<DSMatrixd, DVectord> > &linear_solver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>()
			)
		: B(B), D(D), linear_solver(linear_solver)
		{
			assert(B);
			assert(D);
			assert(linear_solver);

			linear_solver->update(D);
		}

		inline void describe(std::ostream &os) const override
		{
			DVectord t_from = local_values(local_size(*B).get(1), 1);
			DVectord t_to;
			apply(t_from, t_to);

			double t_max = max(t_to);
			double t_min = min(t_to);

			double sum_D = sum(*D);
			double sum_B = sum(*B);

			os << "------------------------------------------\n";
			os << "L2TransferOperator:\n";
			os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
			os << "sum(B) = " << sum_B << ", sum(D) = " << sum_D << std::endl;
			os << "------------------------------------------\n";
		}

		bool write(const Path &path) const override
		{ 
			return utopia::write(path / "B.m", *B) && utopia::write(path / "D.m", *D);
		}

	private:
		std::shared_ptr<DSMatrixd> B;
		std::shared_ptr<DSMatrixd> D;
		std::shared_ptr<LinearSolver<DSMatrixd, DVectord> > linear_solver;
	};

	class PseudoL2TransferOperator : public TransferOperator {
	public:
		inline void apply(const DVectord &from, DVectord &to) const override
		{
			assert(T);
			to = *T * from;
		}

		inline void apply_transpose(const DVectord &from, DVectord &to) const override
		{
			assert(T);
			to = transpose(*T) * from;
		}

		PseudoL2TransferOperator() {}

		inline void init_from_coupling_operator(const DSMatrixd &B)
		{
			T = std::make_shared<DSMatrixd>();
			DVectord d = sum(B, 1);

			{
				ReadAndWrite<DVectord> rw_(d);
				auto r = range(d);
				for(auto k = r.begin(); k != r.end(); ++k) {
					if(approxeq(d.get(k), 0.0, 1e-14)) {
						d.set(k, 1.);
					}
				}
			}

			*T = diag(1./d) * B;
		}

		PseudoL2TransferOperator(const std::shared_ptr<DSMatrixd> &T)
		: T(T)
		{
			assert(T);
		}

		inline void describe(std::ostream &os) const override
		{
			DVectord t = sum(*T, 1);
			double t_max = max(t);
			double t_min = min(t);
			double t_sum = sum(t);

			os << "------------------------------------------\n";
			os << "PseudoL2TransferOperator:\n";
			os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
			os << "sum(T): "  << t_sum << " <= " << size(*T).get(0) << "\n";
			os << "------------------------------------------\n";
		}

		bool write(const Path &path) const override
		{ 
			return utopia::write(path / "T.m", *T);
		}

	private:
		std::shared_ptr<DSMatrixd> T;
	};

	class Interpolator : public TransferOperator {
	public:
		inline void apply(const DVectord &from, DVectord &to) const override
		{
			to = *T * from;
		}

		void apply_transpose(const DVectord &from, DVectord &to) const override
		{
			assert(T);
			to = transpose(*T) * from;
		}

		Interpolator(const std::shared_ptr<DSMatrixd> &T)
		: T(T)
		{
			assert(T);
		}

		void normalize_rows()
		{
			DVectord d = sum(*T, 1);
			ReadAndWrite<DVectord> rw_(d);
			auto r = range(d);
			for(auto k = r.begin(); k != r.end(); ++k) {
				if(approxeq(d.get(k), 0.0, 1e-14)) {
					d.set(k, 1.);
				}
			}

			*T = diag(1./d) * (*T);
		}

		inline void describe(std::ostream &os) const override
		{
			DVectord t = sum(*T, 1);
			double t_max = max(t);
			double t_min = min(t);
			double t_sum = sum(t);

			os << "------------------------------------------\n";
			os << "Interpolator:\n";
			os << "row sum [" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
			os << "sum(T): " << t_sum << " <= " << size(*T).get(0) << "\n";
			os << "------------------------------------------\n";
		}

		bool write(const Path &path) const override
		{ 
			return utopia::write(path / "T.m", *T);
		}

	private:
		std::shared_ptr<DSMatrixd> T;
	};


	enum TransferOperatorType {
		INTERPOLATION = 0,
		L2_PROJECTION = 1,
		PSEUDO_L2_PROJECTION = 2,
		APPROX_L2_PROJECTION = 3
	};
}

#endif //UTOPIA_TRANSFER_ASSEMBLER_HPP
