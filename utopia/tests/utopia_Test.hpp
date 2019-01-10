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
#include "utopia_TaoSolverTest.hpp"
#include "utopia_PetscCudaTest.hpp"
#include "utopia_SelectionTest.hpp"
#include "utopia_UITest.hpp"
#include "utopia_M3ELinSolTest.hpp"
#include "utopia_QPSolverTest.hpp"
#include "utopia_trilinos_KokkosTest.hpp"
#include "utopia_TRTest.hpp"

namespace utopia
{
    inline static void runAllTests()
    {
        runWrapperTest();
        run_autodiff_test();
        runAlgebraTest();
        runUtilitiesTest();
        runPetscTest();
        runBLASTest();
        runMiscTest();
        run_kokkos_test();
        run_trilinos_test();
        run_tao_solver_test();
        run_petsc_cuda_test();
        run_selection_test();
        run_ui_test();


        runGenericSolversTest();
        runPetscNonlinearSolversTest();
        runPetscLinearSolversTest();
        runPetscSlepcSolversTest();
        runQuasiNewtonTest(); 

        runNonlinearMultilevelSolverTest();

        
        run_qp_solver_test();

        //only works for serial
        if(mpi_world_size() == 1) {
            // run_performance_test();
            run_m3e_lin_sol_test();
        }
    }

    inline static void runTests(const std::string& tests)
    {
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin testing]" << std::endl;
        }

        Chrono c;
        c.start();

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
                {
                    runGenericSolversTest();
                    runPetscNonlinearSolversTest();
                    runPetscLinearSolversTest();
                    runPetscSlepcSolversTest();
                }
                else if (token == "solvers_generic")
                   runGenericSolversTest();
                else if (token == "solvers_petsc_nonlinear")
                   runPetscNonlinearSolversTest();
                else if (token == "solvers_petsc_linear")
                   runPetscLinearSolversTest();
                else if (token == "solvers_slepc")
                   runPetscSlepcSolversTest();
                else if(token == "nonlinear_multilevel")
                    runNonlinearMultilevelSolverTest();
                else if(token =="quasi_newton")
                    runQuasiNewtonTest(); 
                // else if (token == "performance")
                    // run_performance_test();
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
                else if (token == "kokkos")
                    run_kokkos_test();
                else if (token == "trilinos")
                    run_trilinos_test();
                else if(token == "tao") {
                    run_tao_solver_test();
                } else if(token == "petsc_cuda") {
                    run_petsc_cuda_test();
                } else if(token == "selection") {
                    run_selection_test();
                } else if(token == "ui") {
                    run_ui_test();
                } else if(token == "m3e") {
                    run_m3e_lin_sol_test();
                } else if(token == "qp") {
                    run_qp_solver_test();
                }else if(token == "uncon_bench") {
                    run_unconstrained_optimization_benchmark();
                }
            }
        }

        if(mpi_world_rank() == 0) {
            std::cout << "[End testing]" << std::endl;
        }

        mpi_world_barrier();
        c.stop();
        if(utopia::Utopia::instance().verbose() && mpi_world_rank() == 0) {
            std::cout << c << std::endl;
        }
    }
}

#endif //UTOPIA_UTOPIA_TEST_HPP
