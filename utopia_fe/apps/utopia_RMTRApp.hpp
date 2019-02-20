#ifndef UTOPIA_RMTR_APP_HPP
#define UTOPIA_RMTR_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "libmesh/parallel_mesh.h"

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	class RMTRApp final : public FEApp {
	public:
		void run(Input &in) override;
		void init(libMesh::Parallel::Communicator &comm) override;

		
		inline static std::string command()
		{
			return "-rmtr";
		}

		class SimulationInput;

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;

		void solve_newton(const SimulationInput &in);
		void solve_rmtr(const SimulationInput &in);
	};
}

#endif //UTOPIA_RMTR_APP_HPP
