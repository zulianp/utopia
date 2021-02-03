
#include <iostream>
#include "utopia.hpp"
#include "utopia_Version.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_impl.hpp"

using VectorType = utopia::PetscVector;
using MatrixType = utopia::PetscMatrix;

#else
#ifdef UTOPIA_WITH_TRILINOS
using VectorType = utopia::TpetraVector;
using MatrixType = utopia::TpetraMatrix;
#else
#ifdef UTOPIA_WITH_BLAS
using VectorType = utopia::BlasVectord;
using MatrixType = utopia::BlasMatrixd;
#else
#error "No backend available"
#endif  // UTOPIA_WITH_BLAS
#endif  // UTOPIA_WITH_TRILINOS
#endif  // UTOPIA_WITH_PETSC

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

    // FIXME use inputs
    auto ls = linear_solver<MatrixType, VectorType>(Solver::automatic());
    if (ls) {
        (*solver)->ptr = ls->clone();
    } else {
        (*solver)->ptr = nullptr;
    }
}

void USolverPrintInfo(USolver solver) {
    using namespace utopia;

    auto solver_ptr = reinterpret_cast<LinearSolver<MatrixType, VectorType> *>(solver->ptr);
    if (solver_ptr) {
        utopia::out() << "USolverPrintInfo: not null\n";
    } else {
        utopia::out() << "USolverPrintInfo: null\n";
    }
    // solver_ptr->describe(std::cout);
}

void USolverDestroy(USolver *solver) {
    using namespace utopia;
    utopia::out() << "USolverDestroy" << std::endl;

    auto solver_ptr = reinterpret_cast<LinearSolver<MatrixType, VectorType> *>((*solver)->ptr);
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
