#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_ILU.hpp"
#include "utopia_ILU_impl.hpp"

#include "utopia_assemble_laplacian_1D.hpp"

using namespace utopia;

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_ILUDecompose.hpp"

void petsc_ilu_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 20;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);
    assemble_laplacian_1D(A, true);

    PetscVector x(vl, 0.0), b(vl, 1.0);
    ILU<PetscMatrix, PetscVector> ls;
    ls.verbose(true);
    ls.atol(1e-6);
    ls.solve(A, b, x);

    if (comm.size() == 1) {
        rename("a", A);
        write("A.m", A);
    }
}

void petsc_ilu_cg_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 20;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);
    assemble_laplacian_1D(A, true);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    auto ilu = std::make_shared<ILU<PetscMatrix, PetscVector>>();
    ilu->max_it(10);

    ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> ls;
    ls.apply_gradient_descent_step(true);

    ls.verbose(true);
    ls.atol(1e-6);

    ls.set_preconditioner(ilu);
    ls.solve(A, b, x);
}

#endif  // UTOPIA_WITH_PETSC

void ilu() {
#ifdef UTOPIA_WITH_PETSC
    UTOPIA_RUN_TEST(petsc_ilu_test);
    UTOPIA_RUN_TEST(petsc_ilu_cg_test);
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(ilu);
