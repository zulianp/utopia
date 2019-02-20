#ifndef UTOPIA_FE_APP_HPP
#define UTOPIA_FE_APP_HPP

#include "utopia_ui.hpp"
#include "utopia_App.hpp"
#include "libmesh/parallel_mesh.h"
#include <string>

namespace utopia {

	class FEApp : public App {
	public:
		virtual ~FEApp() {}
		virtual void init(libMesh::Parallel::Communicator &comm) = 0;
	};
}

#endif //UTOPIA_FE_APP_HPP
