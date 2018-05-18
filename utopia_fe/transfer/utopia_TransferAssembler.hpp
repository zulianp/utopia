#ifndef UTOPIA_TRANSFER_ASSEMBLER_HPP
#define UTOPIA_TRANSFER_ASSEMBLER_HPP 

#include "utopia_LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_libmesh.hpp"

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
			virtual bool assemble(
				LocalAssembler &assembler,
				Local2Global   &local2global,
				SparseMatrix   &B) = 0;
		};

		TransferAssembler();
		~TransferAssembler();

		bool assemble(
			const std::shared_ptr<MeshBase> &from_mesh,
			const std::shared_ptr<DofMap>   &from_dofs,
			const std::shared_ptr<MeshBase> &to_mesh,
			const std::shared_ptr<DofMap>   &to_dofs,
			SparseMatrix &B,
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

	private:
		std::shared_ptr<LocalAssembler> assembler_;
		std::shared_ptr<Local2Global> local2global_;
		std::shared_ptr<Algorithm> algorithm_;
	};

	class TransferOperator {
	public:
		virtual ~TransferOperator() {}
		virtual void apply(const DVectord &from, const DVectord &to) = 0;
	};

	class L2TransferOperator : public TransferOperator {
	public:
		void apply(const DVectord &from, const DVectord &to) override;
	};

	class PseudoL2TransferOperator : public TransferOperator {
	public:
		void apply(const DVectord &from, const DVectord &to) override;
	};

	class Interpolator : public TransferOperator {
	public:
		void apply(const DVectord &from, const DVectord &to) override;
	};
}

#endif //UTOPIA_TRANSFER_ASSEMBLER_HPP
