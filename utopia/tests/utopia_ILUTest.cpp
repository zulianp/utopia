#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_ILU.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

using namespace utopia;

#ifdef UTOPIA_WITH_PETSC

#include "utopia_ILUDecompose.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"
#include "utopia_petsc_ILUDecompose.hpp"

#include "utopia_DILUDecompose_impl.hpp"

void petsc_ilu_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 100;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    assemble_laplacian_1D(A, true);
    // assemble_poisson_problem_1D(1.0, A, b);

    ILU<PetscMatrix, PetscVector> ls;
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

void petsc_dilu_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 1000;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 3, 2);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    assemble_laplacian_1D(A, true);
    // assemble_poisson_problem_1D(1.0, A, b);
    using CrsMatrix_t = CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, 1>;
    using Vector_t = ArrayView<PetscScalar>;

    CrsMatrix_t crs;
    crs_block_matrix(A, crs);

    DILUAlgorithm<CrsMatrix_t, Vector_t> dilu;
    dilu.update(crs);

    // ILU<PetscMatrix, PetscVector> ls;
    // // ls.verbose(true);
    // ls.atol(1e-6);
    // ls.rtol(1e-7);
    // ls.stol(1e-7);
    // ls.max_it(100);

    // ls.solve(A, b, x);

    // PetscVector r = b - A * x;
    // PetscScalar norm_r = norm1(r);
    // // disp(norm_r);

    // utopia_test_assert(norm_r < 1e-5);
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
    ilu->max_it(5);
    // ilu->verbose(true);

    ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> ls;
    ls.apply_gradient_descent_step(true);

    // ls.verbose(true);
    ls.atol(1e-6);
    ls.rtol(1e-7);
    ls.stol(1e-7);

    // if (comm.size() == 1) {
    //     rename("a", A);
    //     write("A.m", A);

    //     rename("b", b);
    //     write("B.m", b);
    // }

    ls.set_preconditioner(ilu);
    ls.update(make_ref(A));
    ls.apply(b, x);

    PetscVector r = b - A * x;
    PetscScalar norm_r = norm1(r);
    utopia_test_assert(norm_r < 1e-5);
}

template <class Matrix>
void assemble_vector2_coupled_laplacian_1D(Matrix &m, const bool bc = false) {
    // n x n matrix with maximum 3 entries x row
    Write<Matrix> w(m);
    Range r = row_range(m);
    auto n = size(m).get(0);

    const double block_d[2 * 2] = {6.0, 2.0, 2.0, 6.0};
    const double block_o[2 * 2] = {-1.0, -1.0, -1.0, -1.0};

    for (SizeType i = r.begin(); i != r.end(); i += 2) {
        if (bc && (i == 0 || i == n - 2)) {
            m.set(i, i, 1.0);
            m.set(i + 1, i + 1, 1.0);
            continue;
        }

        for (SizeType b_i = 0; b_i < 2; ++b_i) {
            for (SizeType b_j = 0; b_j < 2; ++b_j) {
                m.set(i + b_i, i + b_j, block_d[b_i * 2 + b_j]);

                if (i >= 2) {
                    m.set(i + b_i, i - 2 + b_j, block_o[b_i * 2 + b_j]);
                }

                if (i < n - 2) {
                    m.set(i + b_i, i + 2 + b_j, block_o[b_i * 2 + b_j]);
                }
            }
        }
    }
}

void petsc_block_ilu_test() {
    auto comm = PetscCommunicator::get_default();
    PetscInt n = 2000;

    auto vl = layout(comm, n, n * comm.size());
    PetscMatrix A;
    A.sparse(square_matrix_layout(vl), 6, 4);

    PetscVector x(vl, 0.0), b(vl, 1.0);

    assemble_vector2_coupled_laplacian_1D(A, true);
    // assemble_poisson_problem_1D(1.0, A, b);

    InputParameters params;
    params.set("block_size", 2);
    // params.set("print_matrices", true);
    ILU<PetscMatrix, PetscVector> ls;
    ls.read(params);
    // auto ilu = std::make_shared<ILU<PetscMatrix, PetscVector>>();
    // ilu->read(params);
    // ilu->max_it(5);

    // ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE> ls;
    // ls.apply_gradient_descent_step(true);
    // ls.set_preconditioner(ilu);
    ls.verbose(true);
    ls.atol(1e-6);
    ls.rtol(1e-8);
    ls.stol(1e-8);
    ls.max_it(n);
    ls.update(make_ref(A));
    ls.apply(b, x);

    PetscVector r = b - A * x;
    PetscScalar norm_r = norm1(r);

    utopia_test_assert(norm_r < 1e-5);
}

#endif  // UTOPIA_WITH_PETSC

void ilu() {
#ifdef UTOPIA_WITH_PETSC
    UTOPIA_RUN_TEST(petsc_ilu_test);
    UTOPIA_RUN_TEST(petsc_ilu_cg_test);
    UTOPIA_RUN_TEST(petsc_block_ilu_test);
    UTOPIA_RUN_TEST(petsc_dilu_test);
    // UTOPIA_RUN_TEST(petsc_ilu_vi_test);
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(ilu);
