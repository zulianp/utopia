#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_ILU.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

using namespace utopia;

#ifdef UTOPIA_WITH_PETSC

#include "utopia_Agglomerate.hpp"
#include "utopia_ILUDecompose.hpp"
#include "utopia_petsc_ILUDecompose.hpp"

void petsc_amg_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 1000;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    assemble_laplacian_1D(A, true);

    // auto smoother = std::make_shared<GaussSeidel<PetscMatrix, PetscVector>>();
    auto smoother = std::make_shared<ILU<PetscMatrix, PetscVector>>();
    auto coarse_solver = std::make_shared<ILU<PetscMatrix, PetscVector>>();
    // auto coarse_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
    auto agglomerator = std::make_shared<Agglomerate<PetscMatrix>>();

    AlgebraicMultigrid<PetscMatrix, PetscVector> ls(smoother, coarse_solver, agglomerator);
    InputParameters params;
    params.set("n_levels", 5);

    ls.read(params);

    // ls.verbose(true);
    ls.atol(1e-6);
    ls.rtol(1e-7);
    ls.stol(1e-7);
    ls.max_it(100);

    ls.solve(A, b, x);

    PetscVector r = b - A * x;
    PetscScalar norm_r = norm1(r);
    // disp(norm_r);

    utopia_test_assert(norm_r < 1e-5);
}

#endif  // UTOPIA_WITH_PETSC

void amg() {
#ifdef UTOPIA_WITH_PETSC
    UTOPIA_RUN_TEST(petsc_amg_test);
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(amg);
