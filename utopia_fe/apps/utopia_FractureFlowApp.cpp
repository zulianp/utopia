#include "utopia_FractureFlowApp.hpp"
#include "utopia_SteadyFractureFlowSimulation.hpp"
#include "utopia_FractureFlowTransportSimulation.hpp"

namespace utopia {

    void FractureFlowApp::run(Input &in)
    {
        // <flow-type>transient</flow-type>
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
    }
}

