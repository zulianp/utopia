
#include "utopia.hpp"
#include "utopia_AppRunner.hpp"

#include <memory>
#include <iostream>
#include <sstream>
#include <ctime>

using namespace std;

int main(const int argc, char *argv[])
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    {
        //REMOVEME
        for(int i = 1; i < argc; i++) {
            if (argv[i] == std::string("-data_path")) {
                if (++i >= argc)
                    break;
                if (mpi_world_rank() == 0) {
                    std::cout << "data_path: " << argv[i] << std::endl;
                }

                Utopia::instance().set("data_path", argv[i]);
            }
        }

        AppRunner app_runner;
        app_runner.run(argc, argv);
    }

    return Utopia::Finalize();
}
