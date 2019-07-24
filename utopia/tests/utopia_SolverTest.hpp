#ifndef UTOPIA_SOLVER_TEST_HPP
#define UTOPIA_SOLVER_TEST_HPP

namespace utopia
{
    void runGenericSolversTest();
    void runPetscNonlinearSolversTest();
    void runPetscLinearSolversTest();
    void runPetscSlepcSolversTest();
    void runSolversTest();
    void runNonlinearMultilevelSolverTest();
    void runQuasiNewtonTest();
    void runPetscPseudoTransientContinuationTest(); 

    // To be renamed ...
    void runHckTest(); 

}

#endif //UTOPIA_SOLVER_TEST_HPP
