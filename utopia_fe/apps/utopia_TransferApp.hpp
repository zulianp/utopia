#ifndef UTOPIA_TRANSFER_APP
#define UTOPIA_TRANSFER_APP

#include "utopia_App.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"

#include <string>
#include <memory>


namespace utopia {
	class LocalAssembler;
	class Local2Global;

	class TransferApp final : public App {
	public:
		~TransferApp();
		TransferApp();

		void init(libMesh::LibMeshInit &init);
		void run(const std::string &conf_file_path);

		static std::string command() { return "-transfer"; }

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;

		std::shared_ptr<libMesh::DistributedMesh> mesh_master_;
		std::shared_ptr<libMesh::EquationSystems> equation_systems_master_;
		std::shared_ptr<LibMeshFunctionSpace> space_master_;

		std::shared_ptr<libMesh::DistributedMesh> mesh_slave_;
		std::shared_ptr<libMesh::EquationSystems> equation_systems_slave_;
		std::shared_ptr<LibMeshFunctionSpace> space_slave_;

		std::shared_ptr<LocalAssembler> local_assembler_;
		std::shared_ptr<Local2Global> local2global_;
		bool is_interpolation_;

		std::string type;

		std::shared_ptr<TransferOperator> transfer_op_;
		int biorth_basis;


		std::shared_ptr<SymbolicFunction> fun;

	};
}


#endif //UTOPIA_TRANSFER_APP
