#ifndef UTOPIA_FE_APPS_HPP
#define UTOPIA_FE_APPS_HPP

#include "utopia_FEApp.hpp"

#include "libmesh/parallel_mesh.h"
#include <string>


namespace utopia {

	class FEApps {
	public:
		FEApps();
		void run(libMesh::Parallel::Communicator &comm, int argc, char * argv[]);
		void add_app(const std::string &command, std::unique_ptr<FEApp> &&app);

	private:
		std::map<std::string, std::unique_ptr<FEApp>> apps_;
	};
}

#endif //UTOPIA_FE_APPS_HPP
