
#include "utopia.hpp"
#include "utopia_Test.hpp"
#include <memory>
#include <iostream>
#include <sstream>
#include <ctime>

using namespace std;

int main(const int argc, char *argv[])
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    std::string tests = "all";

    for (size_t i = 1; i < argc; i++) {
        if (argv[i] == std::string("-data_path")) {
            if (++i >= argc)
                break;
            if (mpi_world_rank() == 0) {
                std::cout << "data_path: " << argv[i] << std::endl;
            }
            Utopia::Instance().set("data_path", argv[i]);
        } else if (argv[i] == std::string("-test")) {
            if (++i >= argc)
                break;
            tests = argv[i];
        }
    }

    // Utopia::Instance().set("default_tollerance", "1e-15");
    // Utopia::Instance().set("n_threads", "2");

    runTests(tests);
    return Utopia::Finalize();
}
