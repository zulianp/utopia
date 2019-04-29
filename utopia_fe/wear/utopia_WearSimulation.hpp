#ifndef UTOPIA_FE_WEAR_SIMULATION_HPP
#define UTOPIA_FE_WEAR_SIMULATION_HPP

#include <string>

namespace libMesh {
    class LibMeshInit;
}

namespace utopia {
    class WearSimulation {
    public:
        class SimulationInput;
        class SimulationOutput;

        void run(libMesh::LibMeshInit &init, const std::string &conf_file_path);

        WearSimulation();
        ~WearSimulation();
    };
}

#endif //UTOPIA_FE_WEAR_SIMULATION_HPP
