#include "utopia_FEApps.hpp"
#include "utopia_WearSimulation.hpp"
#include "utopia_TransferApp.hpp"
#include "utopia_FractureFlowApp.hpp"
#include "utopia_RMTRApp.hpp"
#include "utopia_Grid2MeshTransferApp.hpp"
#include "utopia_ContactApp.hpp"
#include "utopia_EikonalApp.hpp"
#include "utopia_LeastSquaresHelmholtzApp.hpp"
#include "utopia_PoissonApp.hpp"
#include "utopia_SmoothApp.hpp"
#include "utopia_make_unique.hpp"

#include <iostream>

namespace utopia {

    void FEApps::print_usage(std::ostream &os) const
    {
        os << "------------------------------------------------------\n";
        os << "available apps:\n";

        for(const auto &t : apps_) {
            os << "\t" << t.first << "\n";
        }

        os << std::endl;
        os << "------------------------------------------------------\n";
    }


    void FEApps::add_app(const std::string &command, std::unique_ptr<FEApp> &&app)
    {
        apps_[command] = std::move(app);
    }

    void FEApps::run(libMesh::Parallel::Communicator &comm, int argc, char * argv[])
    {
        for(int i = 1; i < argc; ++i) {
            const int ip1 = i+1;

            if(argv[i] == std::string("-test")) {
                i = ip1;
                continue;
            }

            auto it = apps_.find(argv[i]);
            if(it == apps_.end()) continue;

            if(ip1 >= argc) {
                std::cerr << "[Error] expected input file (xml, json)" << std::endl;
                return;
            }

            auto in = open_istream(argv[ip1]);

            if(!in) {
                std::cerr << "[Error] no valid input file found at path: " << argv[ip1] << std::endl;
                return;
            }

            it->second->init(comm);
            it->second->run(*in);
        }
    }

    FEApps::FEApps()
    {
        add_app(TransferApp::command(),  		 	 utopia::make_unique<TransferApp>());
        add_app(FractureFlowApp::command(),  	 	 utopia::make_unique<FractureFlowApp>());
        add_app(RMTRApp::command(),  			 	 utopia::make_unique<RMTRApp>());
        add_app(ContactApp::command(),  		 	 utopia::make_unique<ContactApp>());
        add_app(Grid2MeshTransferApp::command(), 	 utopia::make_unique<Grid2MeshTransferApp>());
        add_app(EikonalApp::command(), 			 	 utopia::make_unique<EikonalApp>());
        add_app(LeastSquaresHelmholtzApp::command(), utopia::make_unique<LeastSquaresHelmholtzApp>());
        add_app(PoissonApp::command(), utopia::make_unique<PoissonApp>());
        add_app(SmoothApp::command(), utopia::make_unique<SmoothApp>());
    }
}
