#ifndef UTOPIA_FE_WEAR_SIMULATION_HPP
#define UTOPIA_FE_WEAR_SIMULATION_HPP

#include <string>

namespace libMesh {
	class LibMeshInit;
}

namespace utopia {
	class WearSimulation {
	public:
		class Impl;
		void run(libMesh::LibMeshInit &init, const std::string &xml_file_path);

		WearSimulation();
		~WearSimulation();
	};
}

#endif //UTOPIA_FE_WEAR_SIMULATION_HPP
