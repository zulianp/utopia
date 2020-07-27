
#include <iostream>
#include "utopia_Instance.hpp"
#include "utopia_Version.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"

extern "C" {

#include "utopia_c_solve.h"

void UtopiaPrintVersion() { utopia::out() << UTOPIA_VERSION << std::endl; }

void UtopiaInitialize(int argc, char *argv[]) {
    using namespace utopia;
    utopia::out() << "UtopiaInitialize" << std::endl;
    Utopia::Init(argc, argv);
}

int UtopiaFinalize() {
    utopia::out() << "UtopiaFinalize" << std::endl;

    using namespace utopia;
    return Utopia::Finalize();
}

void USolverCreate(USolver *solver, USolverType type, UPreconditionerType prec, UPackage package) {
    // FIXME
    UTOPIA_UNUSED(prec);
    UTOPIA_UNUSED(package);

    using namespace utopia;
    utopia::out() << "USolverCreate: (USolverType=" << type << ")" << std::endl;
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
        utopia::out() << "KSP" << std::endl;
        ksp->describe(std::cout);
    } else {
        auto fact = dynamic_cast<Factorization<PetscMatrix, PetscVector> *>(solver_ptr);
        fact->describe(std::cout);
    }
}

void USolverDestroy(USolver *solver) {
    using namespace utopia;
    utopia::out() << "USolverDestroy" << std::endl;

    auto solver_ptr = reinterpret_cast<LinearSolver<PetscMatrix, PetscVector> *>((*solver)->ptr);
    delete solver_ptr;
    delete *solver;
    *solver = nullptr;
}

void USolverSolve(USolver solver, UMat A, UVec b, UVec x) {
    // FIXME
    UTOPIA_UNUSED(solver);
    UTOPIA_UNUSED(A);
    UTOPIA_UNUSED(b);
    UTOPIA_UNUSED(x);

    utopia::out() << "USolverSolve" << std::endl;
}
}
