#ifndef UTOPIA_TRANSFER_APP
#define UTOPIA_TRANSFER_APP

#include "utopia_App.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_libmesh.hpp"

#include <string>
#include <memory>


namespace utopia {
	class LocalAssembler;
	class Local2Global;

	class TransferApp final : public App {
	public:
		class InputSpace;

		~TransferApp();
		TransferApp();

		void init(libMesh::LibMeshInit &init);
		void run(const std::string &conf_file_path);

		static std::string command() { return "-transfer"; }

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;

		std::shared_ptr<LocalAssembler> local_assembler_;
		std::shared_ptr<Local2Global> local2global_;
		bool is_interpolation_;

		std::string type;

		std::shared_ptr<TransferOperator> transfer_op_;
		int biorth_basis;
		int assemble_mass_mat_;


#ifdef WITH_TINY_EXPR
		std::shared_ptr<SymbolicFunction> fun;
#else
		std::shared_ptr<ConstantCoefficient<double, 0>> fun;
#endif //WITH_TINY_EXPR

		bool fun_is_constant;
		bool write_operators_to_disk;

	};
}


#endif //UTOPIA_TRANSFER_APP
