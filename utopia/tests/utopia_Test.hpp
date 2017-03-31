#ifndef UTOPIA_UTOPIA_TEST_HPP
#define UTOPIA_UTOPIA_TEST_HPP

#include <sstream>
#include "utopia_MemoryPool.hpp"
#include "utopia_SpecTest.hpp"
#include "utopia_WrapperTest.hpp"
#include "test_problems/utopia_TestProblems.hpp"
#include "utopia_AutoDiffTest.hpp"
#include "utopia_SolverTest.hpp"
#include "utopia_PerformanceTest.hpp"
#include "utopia_AlgebraTest.hpp"
#include "utopia_UtilitiesTest.hpp"
#include "utopia_PETScTest.hpp"
#include "utopia_BLASTest.hpp"
#include "utopia_MiscTest.hpp"

namespace utopia
{
    inline static void runAllTests()
    {
        runWrapperTest();
        runSpecTest();
        run_autodiff_test();
        runSolversTest();
        runAlgebraTest();
        runUtilitiesTest();
        runPETScTest();
        runBLASTest();
        runMiscTest();


        //only works for serial
        if(mpi_world_size() == 1) {
            run_performance_test();
        }
    }

    inline static void runTests(const std::string& tests)
    {
        if (tests == "all") {
            runAllTests();
            return;
        }
        std::istringstream iss(tests);
        std::string token;
        while (std::getline(iss, token, ',')) {
            MEMPOOL().fullGC();
            if (token == "wrapper")
                runWrapperTest();
            else if (token == "spec")
                runSpecTest();
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
                runPETScTest();
            else if (token == "blas")
                runBLASTest();
            else if (token == "misc")
                runMiscTest();
        }
    }
}

#endif //UTOPIA_UTOPIA_TEST_HPP
