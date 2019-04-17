#include <iostream>

#include "utopia_libmesh.hpp"
#include "par_moonolith.hpp"
#include "utopia.hpp"

#include "utopia_Socket.hpp"

#include "utopia_FEApps.hpp"
#include "utopia_FETests.hpp"

#include <functional>
#include <iostream>

using namespace utopia;
using namespace std;
using namespace libMesh;

int main(int argc, char *argv[])
{
    //For debugginh with ddt
    MPI_Init(&argc, &argv);
    PETSC_COMM_WORLD = MPI_COMM_WORLD;

    Utopia::Init(argc, argv);
    MOONOLITH_PROFILING_BEGIN();

    {
        LibMeshInit init(argc, argv, PETSC_COMM_WORLD);

        FEApps apps;
        FETests tests;

        for(int i = 1; i < argc; ++i) {
            const int ip1 = i+1;

            if(argv[i] == std::string("-verbose")) {
                Utopia::instance().set("verbose", "true");
                continue;
            } else if(argv[i] == std::string("-h")) {
                std::cout << "--------------------------------------------" << std::endl;
                std::cout << "--------------------------------------------" << std::endl;
                apps.print_usage(std::cout);
                tests.print_usage(std::cout);


                std::cout << "--------------------------------------------" << std::endl;
                std::cout << "--------------------------------------------" << std::endl;
            } else if(argv[i] == std::string("-output_path")) {
                utopia::Utopia::instance().set("output_path", argv[ip1]);
                std::cout << "setting output_path to: " << argv[ip1] << std::endl;
            } else if(argv[i] == std::string("-data_path")) {
                utopia::Utopia::instance().set("data_path", argv[ip1]);
                std::cout << "setting data_path to: " << argv[ip1] << std::endl;
            }
        }


        apps.run(init.comm(), argc, argv);
        tests.run(init.comm(), argc, argv);
    }

    MOONOLITH_PROFILING_END();
    return Utopia::Finalize();
}


