
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
        AppRunner app_runner;
        std::vector<std::string> apps;
        for (int i = 1; i < argc; i++) {
            if (argv[i] == std::string("-data_path")) {
                if (++i >= argc)
                    break;
                if (mpi_world_rank() == 0) {
                    std::cout << "data_path: " << argv[i] << std::endl;
                }

                Utopia::instance().set("data_path", argv[i]);
            } else if (argv[i] == std::string("-app")) {
                if (++i >= argc)
                    break;
                apps.push_back(argv[i]);
            } else if(argv[i] == std::string("-list")) {
                app_runner.describe();
            }
        }

        app_runner.verbose(Utopia::instance().verbose());
        app_runner.run(apps);
    }

    return Utopia::Finalize();
}
