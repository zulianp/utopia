#ifndef UTOPIA_GRID_2_MESH_TRANSFER_APP
#define UTOPIA_GRID_2_MESH_TRANSFER_APP

#include "utopia_FEApp.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_libmesh.hpp"

#include <string>
#include <memory>


namespace utopia {
	class LocalAssembler;
	class Local2Global;

	class Grid2MeshTransferApp final : public FEApp {
	public:
		class InputSpace;
		class InputGrid;

		~Grid2MeshTransferApp();
		Grid2MeshTransferApp();

		void run(Input &in) override;
		void init(libMesh::Parallel::Communicator &comm) override;


		static std::string command() { return "-g2m_transfer"; }

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;
	};
}


#endif //UTOPIA_GRID_2_MESH_TRANSFER_APP
