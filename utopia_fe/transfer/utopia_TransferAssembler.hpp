#ifndef UTOPIA_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_TRANSFER_ASSEMBLER_HPP

#include "utopia_libmesh.hpp"
#include "utopia_LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia.hpp"
#include "utopia_Path.hpp"

#include <memory>

namespace utopia {

	class TransferOptions {
	public:
		TransferOptions()
		: from_var_num(0), to_var_num(0), n_var(1), to_trace_space(false), tags({})
		{}

		int from_var_num;
		int to_var_num;
		int n_var;
		bool to_trace_space;
		std::vector< std::pair<int, int> > tags;
	};

	class TransferAssembler final {
	public:
		using FunctionSpace = utopia::LibMeshFunctionSpace;
		using SparseMatrix  = utopia::USparseMatrix;
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
		virtual void apply(const UVector &from, UVector &to) const = 0;
		virtual void apply_transpose(const UVector &from, UVector &to) const = 0;
		virtual void describe(std::ostream &) const {}
		virtual bool write(const Path &) const { return false; }
	};


	/**
	 * @brief constructed as (D^-1 * B) * ( . )
	 */
	class L2TransferOperator final : public TransferOperator {
	public:
		inline void apply(const UVector &from, UVector &to) const override
		{
			UVector B_from = *B * from;

			if(empty(to)) {
				to = B_from;
			}

			linear_solver->apply(B_from, to);
		}

		void fix_mass_matrix_operator()
		{
			UVector d;

			Size s = local_size(*D);
			d = local_values(s.get(0), 1.);

			{
				Write<UVector> w_d(d);

				each_read(*D, [&d](const SizeType i, const SizeType, const double) {
					d.set(i, 0.);
				});
			}

			(*D) += USparseMatrix(diag(d));
		}

		///@brief assumes that D is symmetric
		void apply_transpose(const UVector &from, UVector &to) const override
		{
			UVector D_inv_from = local_zeros(local_size(*D).get(0));
			linear_solver->apply(from, D_inv_from);
			to = transpose(*B) * D_inv_from;
		}

		inline L2TransferOperator(
			const std::shared_ptr<USparseMatrix> &B,
			const std::shared_ptr<USparseMatrix> &D,
			const std::shared_ptr<LinearSolver<USparseMatrix, UVector> > &linear_solver = std::make_shared<BiCGStab<USparseMatrix, UVector>>()
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
			UVector t_from = local_values(local_size(*B).get(1), 1);
			UVector t_to;
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
		std::shared_ptr<USparseMatrix> B;
		std::shared_ptr<USparseMatrix> D;
		std::shared_ptr<LinearSolver<USparseMatrix, UVector> > linear_solver;
	};

	class PseudoL2TransferOperator final : public TransferOperator {
	public:
		inline void apply(const UVector &from, UVector &to) const override
		{
			assert(T);
			to = *T * from;
		}

		inline void apply_transpose(const UVector &from, UVector &to) const override
		{
			assert(T);
			to = transpose(*T) * from;
		}

		PseudoL2TransferOperator() {}

		inline void init_from_coupling_operator(const USparseMatrix &B)
		{
			T = std::make_shared<USparseMatrix>();
			UVector d = sum(B, 1);

			{
				ReadAndWrite<UVector> rw_(d);
				auto r = range(d);
				for(auto k = r.begin(); k != r.end(); ++k) {
					if(approxeq(d.get(k), 0.0, 1e-14)) {
						d.set(k, 1.);
					}
				}
			}

			*T = diag(1./d) * B;
		}

		PseudoL2TransferOperator(const std::shared_ptr<USparseMatrix> &T)
		: T(T)
		{
			assert(T);
		}

		inline void describe(std::ostream &os) const override
		{
			UVector t = sum(*T, 1);
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

		inline std::shared_ptr<USparseMatrix> matrix() {
			return T;
		}

	private:
		std::shared_ptr<USparseMatrix> T;
	};

	// class PermutedOperator final : public TransferOperator {
	// public:
	// 	PermutedOperator(
	// 		const std::shared_ptr<TransferOperator> &op,
	// 		const std::shared_ptr<USparseMatrix> &from_permutation
	// 		const std::shared_ptr<USparseMatrix> &to_permutation)
	// 	: op_(op), from_permutation_(from_permutation), to_permutation_(to_permutation)
	// 	{}

	// 	inline void apply(const UVector &from, UVector &to) const
	// 	{
	// 		if(from_permutation_) {
	// 			from_buffer_ = *from_permutation_ * from;
	// 		} else {
	// 			from_buffer_ = from;
	// 		}
			
	// 		op_->apply(from_buffer_, to_buffer_);

	// 		if(to_permutation_) {
	// 			to = transpose(*to_permutation_) * to_buffer_;
	// 		} else {
	// 			to = to_buffer_;
	// 		}
	// 	}

	// 	inline void apply_transpose(const UVector &to, UVector &from) const
	// 	{
	// 		if(to_permutation_) {
	// 			to_buffer_ = *to_permutation_ * to;
	// 		} else {
	// 			to_buffer_ = to;
	// 		}
			
	// 		op_->apply_transpose(to_buffer_, from_buffer);

	// 		if(from_permutation_) {
	// 			from = transpose(*from_permutation_) * from_buffer;
	// 		} else {
	// 			from = from_buffer_;
	// 		}
	// 	}


	// 	std::shared_ptr<TransferOperator> op_;
	// 	std::shared_ptr<USparseMatrix> from_permutation_;
	// 	std::shared_ptr<USparseMatrix> to_permutation_;

	// 	UVector from_buffer_, to_buffer_;

	// };

	class BidirectionalOperator final : public TransferOperator {
	public:
		BidirectionalOperator(
			const std::shared_ptr<TransferOperator> &forward,
			const std::shared_ptr<TransferOperator> &backward
		) : forward_(forward), backward_(backward)
		{}

		inline void apply(const UVector &from, UVector &to) const
		{
			forward_->apply(from, to);
		}

		inline void apply_transpose(const UVector &from, UVector &to) const
		{
			backward_->apply(from, to);
		}

		inline void describe(std::ostream &os) const 
		{
			forward_->describe(os);
			backward_->describe(os);
		}

		inline bool write(const Path &) const { return false; }

	private:
		std::shared_ptr<TransferOperator> forward_;
		std::shared_ptr<TransferOperator> backward_;
	};

	class Interpolator final : public TransferOperator {
	public:
		inline void apply(const UVector &from, UVector &to) const override
		{
			to = *T * from;
		}

		void apply_transpose(const UVector &from, UVector &to) const override
		{
			assert(T);
			to = transpose(*T) * from;
		}

		Interpolator(const std::shared_ptr<USparseMatrix> &T)
		: T(T)
		{
			assert(T);
		}

		void normalize_rows()
		{
			UVector d = sum(*T, 1);
			ReadAndWrite<UVector> rw_(d);
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
			UVector t = sum(*T, 1);
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
		std::shared_ptr<USparseMatrix> T;
	};


	enum TransferOperatorType {
		INTERPOLATION = 0,
		L2_PROJECTION = 1,
		PSEUDO_L2_PROJECTION = 2,
		APPROX_L2_PROJECTION = 3,
		BIDIRECTIONAL_L2_PROJECTION = 4,
		BIDIRECTIONAL_PSEUDO_L2_PROJECTION = 5,
	};
}

#endif //UTOPIA_TRANSFER_ASSEMBLER_HPP
