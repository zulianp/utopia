#include "utopia_FractureFlowApp.hpp"
#include "utopia_SteadyFractureFlowSimulation.hpp"
#include "utopia_FractureFlowTransportSimulation.hpp"

namespace utopia {
    
    void FractureFlowApp::init(libMesh::LibMeshInit &init)
    {
        comm_ = make_ref(init.comm());
    }

    void FractureFlowApp::run(const std::string &conf_file_path)
    {
        auto is_ptr = open_istream(conf_file_path);

        // <flow-type>transient</flow-type>
        std::string flow_type;

        is_ptr->get("flow-type", flow_type);

        if(flow_type.empty() || flow_type == "steady") {
            SteadyFractureFlowSimulation sim(*comm_);
            sim.read(*is_ptr);
            sim.run();
        } else {
            FractureFlowTransportSimulation sim(*comm_);
            sim.read(*is_ptr);
            sim.run();
        }
    }
}

