
#include <iostream>
#include "utopia_Instance.hpp"
#include "utopia_Version.hpp"
#include "utopia_petsc.hpp"

extern "C" {

#include "utopia_c_solve.h"

void UtopiaPrintVersion() { std::cout << UTOPIA_VERSION << std::endl; }

void UtopiaInitialize(int argc, char *argv[]) {
    using namespace utopia;
    std::cout << "UtopiaInitialize" << std::endl;
    Utopia::Init(argc, argv);
}

int UtopiaFinalize() {
    std::cout << "UtopiaFinalize" << std::endl;

    using namespace utopia;
    return Utopia::Finalize();
}

void USolverCreate(USolver *solver, USolverType type, UPreconditionerType prec, UPackage package) {
    using namespace utopia;
    std::cout << "USolverCreate: (USolverType=" << type << ")" << std::endl;
    *solver = new USolverImpl;

    static const std::string ksp = "ksp";
    static const std::string fact = "fact";

    if (fact == type) {
        (*solver)->ptr = new Factorization<PetscMatrix, PetscVector>();
    } else {
        (*solver)->ptr = new KSPSolver<PetscMatrix, PetscVector>();
    }
}

void USolverPrintInfo(USolver solver) {
    using namespace utopia;

    auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector> *>(solver->ptr);
    auto ksp = dynamic_cast<KSPSolver<PetscMatrix, PetscVector> *>(solver_ptr);

    if (ksp) {
        std::cout << "KSP" << std::endl;
        ksp->describe(std::cout);
    } else {
        auto fact = dynamic_cast<Factorization<PetscMatrix, PetscVector> *>(solver_ptr);
        fact->describe(std::cout);
    }
}

void USolverDestroy(USolver *solver) {
    using namespace utopia;
    std::cout << "USolverDestroy" << std::endl;

    auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector> *>((*solver)->ptr);
    delete solver_ptr;
    delete *solver;
    *solver = nullptr;
}

void USolverSolve(USolver solver, UMat A, UVec b, UVec x) { std::cout << "(void" << std::endl; }
}
