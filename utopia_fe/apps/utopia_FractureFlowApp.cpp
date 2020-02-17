#include "utopia_FractureFlowApp.hpp"
#include "utopia_SteadyFractureFlowSimulation.hpp"
#include "utopia_FractureFlowTransportSimulation.hpp"
#include "utopia_FracturedPourousMedia.hpp"

namespace utopia {

    void FractureFlowApp::run(Input &in)
    {
        // <flow-type>transient</flow-type>
        int version = 1;
        in.get("version", version);

        utopia::Utopia::instance().read(in);

        std::cout << "disable-adaptivity: " << utopia::Utopia::instance().get("disable-adaptivity")  << std::endl;

        if(version == 1) {
            std::string flow_type;
            in.get("flow-type", flow_type);
            if(flow_type.empty() || flow_type == "steady") {
                SteadyFractureFlowSimulation sim(comm());
                sim.read(in);
                sim.run();
            } else {
                FractureFlowTransportSimulation sim(comm());
                sim.read(in);
                sim.run();
            }
        } else if(version == 2) {

            std::string flow_type;
            in.get("flow-type", flow_type);

            if(flow_type == "dfn-only") {
                DiscreteFractureNetwork<USparseMatrix, UVector> dfn(comm());
                in.get("fracture-network", dfn);
                dfn.compute_flow();
            } else {
                FracturedPourousMedia<USparseMatrix, UVector> fpm(comm());
                fpm.read(in);
                fpm.compute_flow();
            }
        }
    }
}

