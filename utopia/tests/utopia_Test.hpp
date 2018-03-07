#ifndef UTOPIA_UTOPIA_TEST_HPP
#define UTOPIA_UTOPIA_TEST_HPP

#include <sstream>
#include "utopia_WrapperTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_AutoDiffTest.hpp"
#include "utopia_SolverTest.hpp"
#include "utopia_PerformanceTest.hpp"
#include "utopia_AlgebraTest.hpp"
#include "utopia_UtilitiesTest.hpp"
#include "utopia_PetscTest.hpp"
#include "utopia_BLASTest.hpp"
#include "utopia_MiscTest.hpp"
#include "utopia_TrilinosTest.hpp"

namespace utopia
{
    inline static void runAllTests()
    {
        runWrapperTest();
        run_autodiff_test();
        runSolversTest();
        runAlgebraTest();
        runUtilitiesTest();
        runPetscTest();
        runBLASTest();
        runMiscTest();
        run_trilinos_test();


        //only works for serial
        if(mpi_world_size() == 1) {
            run_performance_test();
        }
    }

    inline static void runTests(const std::string& tests)
    {
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin testing]" << std::endl;
        }

        if (tests == "all") {
            runAllTests();
        } else {
            std::istringstream iss(tests);
            std::string token;
            while (std::getline(iss, token, ',')) {
                if (token == "wrapper")
                    runWrapperTest();
                else if (token == "autodiff")
                    run_autodiff_test();
                else if (token == "solvers")
                   runSolversTest();
                else if (token == "performance")
                    run_performance_test();
                else if (token == "algebra")
                    runAlgebraTest();
                else if (token == "utilities")
                    runUtilitiesTest();
                else if (token == "petsc")
                    runPetscTest();
                else if (token == "blas")
                    runBLASTest();
                else if (token == "misc")
                    runMiscTest();
                else if (token == "trilinos")
                    run_trilinos_test();
            }
        }

        if(mpi_world_rank() == 0) {
            std::cout << "[End testing]" << std::endl;
        }
    }
}

#endif //UTOPIA_UTOPIA_TEST_HPP
