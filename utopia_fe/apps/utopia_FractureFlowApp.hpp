#ifndef UTOPIA_FRACTURE_FLOW_APP_HPP
#define UTOPIA_FRACTURE_FLOW_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "libmesh/parallel_mesh.h"

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	class FractureFlowApp final : public FEApp {
	public:
		void run(Input &in) override;
		void init(libMesh::Parallel::Communicator &comm) override;

		
		inline static std::string command()
		{
			return "-fracflow";
		}

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;
	};
}

#endif //UTOPIA_FRACTURE_FLOW_APP_HPP
