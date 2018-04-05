
#include "utopia.hpp"
#include "utopia_Test.hpp"
#include "utopia_benchmarks.hpp"

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
    bool run_tests = true;

    for (int i = 1; i < argc; i++) {
        if (argv[i] == std::string("-data_path")) {
            if (++i >= argc)
                break;
            if (mpi_world_rank() == 0) {
                std::cout << "data_path: " << argv[i] << std::endl;
            }
            Utopia::instance().set("data_path", argv[i]);
        } else if (argv[i] == std::string("-test")) {
            if (++i >= argc)
                break;
            tests = argv[i];
        } else if(argv[i] == std::string("-skip-tests")) {
            run_tests = false;
        } else if(argv[i] == std::string("-verbose")) {
            Utopia::instance().set("verbose", "true");
        } else if(argv[i] == std::string("-performance_test_verbose")) {
            Utopia::instance().set("performance_test_verbose", "true");
        } else if(argv[i] == std::string("-bench")) {
            run_benchmarks();
        }
    }

    // Utopia::instance().set("default_tollerance", "1e-15");
    // Utopia::instance().set("n_threads", "2");
   
    if(run_tests) { 
        runTests(tests);
    }
    
    return Utopia::Finalize();
}
