#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_ILU.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

using namespace utopia;

#ifdef UTOPIA_WITH_PETSC

#include "utopia_petsc_ILUDecompose.hpp"

void petsc_ilu_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 1000;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    assemble_laplacian_1D(A, true);
    // assemble_poisson_problem_1D(1.0, A, b);

    ILU<PetscMatrix, PetscVector> ls;
    ls.verbose(true);
    ls.atol(1e-6);
    ls.rtol(1e-6);
    ls.stol(1e-6);
    ls.solve(A, b, x);
}

void petsc_ilu_cg_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 1e2;
    PetscInt n_global = n * comm.size();

    auto vl = layout(comm, n, n_global);
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);

    PetscVector x(vl, 0.0), b(vl, 0.0);
    assemble_poisson_problem_1D(1.0, A, b, false);

    auto ilu = std::make_shared<ILU<PetscMatrix, PetscVector>>();
    ilu->max_it(10);

    ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> ls;
    ls.apply_gradient_descent_step(true);

    ls.verbose(true);
    ls.atol(1e-6);
    ls.rtol(1e-6);

    // if (comm.size() == 1) {
    //     rename("a", A);
    //     write("A.m", A);

    //     rename("b", b);
    //     write("B.m", b);
    // }

    ls.set_preconditioner(ilu);
    ls.update(make_ref(A));
    ls.apply(b, x);
}

// void petsc_ilu_vi_test() {
//     auto comm = PetscCommunicator::get_default();
//     PetscInt n = 20;
//     PetscInt n_global = n * comm.size();

//     auto vl = layout(comm, n, n_global);
//     PetscMatrix A;
//     A.sparse(square_matrix_layout(vl), 3, 2);

//     // auto h = 1. / (n_global - 1);
//     auto h = 1.;
//     assemble_laplacian_1D_with_scaling(A, h, true);

//     PetscVector x(vl, 0.0), b(vl, 1.0 * h);
//     PetscVector c(vl, 0.0), r(vl, 0.0);
//     PetscVector lb(vl, -1), ub(vl, 40);

//     {
//         Write<PetscVector> w(b);
//         auto r = b.range();

//         if (r.inside(0)) {
//             b.set(0, 1.0);
//         }

//         if (r.inside(n_global - 1)) {
//             b.set(n_global - 1, 1.0);
//         }
//     }

//     PetscMatrix ilu;
//     ILUDecompose<PetscMatrix>::decompose(A, ilu, true);

//     PetscVector dlb, dub;
//     for (SizeType i = 0; i < 100; ++i) {
//         dlb = lb - x;
//         dub = ub - x;

//         r = A * x;
//         r = b - r;

//         c.set(0.0);
//         ILUDecompose<PetscMatrix>::apply_vi(ilu, dlb, dub, r, c);
//         // ILUDecompose<PetscMatrix>::apply(ilu, r, c);

//         x += c;

//         PetscScalar norm_c = norm2(c);
//         PetscScalar norm_r = norm2(r);
//         utopia::out() << i << ") norm_r: " << norm_r << ", norm_c: " << norm_c << "\n";

//         if (norm_c < 1e-8) break;
//     }

//     rename("x", x);
//     write("X.m", x);

//     rename("a", A);
//     write("A.m", A);

//     rename("b", b);
//     write("B.m", b);

//     rename("m", ilu);
//     write("M.m", ilu);

//     // disp(x);
//     // auto ilu = std::make_shared<ILU<PetscMatrix, PetscVector>>();
//     // ilu->max_it(10);

//     // ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> ls;
//     // ls.apply_gradient_descent_step(true);

//     // // ls.verbose(true);
//     // ls.atol(1e-6);
//     // ls.rtol(1e-6);

//     // ls.set_preconditioner(ilu);
//     // ls.update(make_ref(A));
//     // ls.apply(b, x);
// }

#endif  // UTOPIA_WITH_PETSC

void ilu() {
#ifdef UTOPIA_WITH_PETSC
    UTOPIA_RUN_TEST(petsc_ilu_test);
    UTOPIA_RUN_TEST(petsc_ilu_cg_test);
    // UTOPIA_RUN_TEST(petsc_ilu_vi_test);
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(ilu);
