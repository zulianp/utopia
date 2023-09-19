#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_petsc_PatchSmoother.hpp"
#include "utopia_petsc_RASPatchSmoother.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include "test_problems/utopia_QPSolverTestProblem.hpp"

using namespace utopia;

template <class Matrix, class Vector>
class PatchSmootherTest {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using Comm = typename Traits::Communicator;

#ifdef UTOPIA_WITH_BLAS
    using PatchMatrix = utopia::BlasMatrixd;
    using PatchVector = utopia::BlasVectord;
#else
    using PatchMatrix = Matrix;
    using PatchVector = Vector;
#endif  // UTOPIA_WITH_BLAS

    SizeType n_dofs{40};
    bool verbose{false};

    void run() {
        UTOPIA_RUN_TEST(test_qp_patch_smoother_RAS);

        UTOPIA_RUN_TEST(test_qp_patch_smoother);
        // UTOPIA_RUN_TEST(test_qp_patch_smoother_ssn);

        // UTOPIA_RUN_TEST(test_linear_patch_smoother);
    }

    void test_linear_patch_smoother() {
        auto solver = std::make_shared<ProjectedGaussSeidel<PatchMatrix, PatchVector>>();

        PatchSmoother<Matrix, PatchMatrix> patch_smoother;
        patch_smoother.set_patch_solver(solver);

        test_1D("patch_smoother", patch_smoother);
    }

    void test_qp_patch_smoother() {
        if(Comm::get_default().size() > 1) return;

        auto solver = std::make_shared<ProjectedGaussSeidel<PatchMatrix, PatchVector>>();
        PatchSmoother<Matrix, PatchMatrix> patch_smoother;
        patch_smoother.set_patch_solver(solver);

        QPSolverTestProblem<Matrix, Vector>::run(n_dofs, verbose, patch_smoother);
    }

    void test_qp_patch_smoother_RAS() {
        if(Comm::get_default().size() > 1) return;

        auto solver = std::make_shared<ProjectedGaussSeidel<PatchMatrix, PatchVector>>();
        RASPatchSmoother<Matrix, PatchMatrix> patch_smoother;
        patch_smoother.set_patch_solver(solver);
        InputParameters params;
        params.set("overlap", 4);
        patch_smoother.read(params);

        // std::string path = Utopia::instance().get("data_path");

        // Matrix K;
        // read(path + "/RHS_10x10x10_hexa_3D", rhs);
        // read(path + "/K_hexa_10x10x10_3D", K);
        // read(path + "/M_hexa_10x10x10_3D", M);

        // patch_smoother.update(make_ref(K));

        QPSolverTestProblem<Matrix, Vector>::run(n_dofs, verbose, patch_smoother);
    }

    void test_qp_patch_smoother_ssn() {
        if(Comm::get_default().size() > 1) return;
        
        auto solver =
            std::make_shared<SemismoothNewton<Matrix, Vector>>(std::make_shared<Factorization<Matrix, Vector>>());

        // solver->verbose(true);

        PatchSmoother<Matrix> patch_smoother;
        patch_smoother.set_patch_solver(solver);

        QPSolverTestProblem<Matrix, Vector>::run(n_dofs, verbose, patch_smoother);
    }

    template <class PatchSmoother>
    void test_1D(const std::string &name, PatchSmoother &patch_smoother) {
        UTOPIA_UNUSED(name);

        auto comm = Comm::get_default();

        auto vl = layout(comm, 10, 10 * comm.size());
        Matrix A;
        A.sparse(square_matrix_layout(vl), 3, 2);

        Vector x(vl, 0.0), b(vl, 0.01);

        assemble_laplacian_1D(A, true);

        auto rr = b.range();

        {
            Write<Vector> w(b);

            if (rr.inside(0)) {
                b.set(0, 0);
            }

            if (rr.inside(b.size() - 1)) {
                b.set(b.size() - 1, 0);
            }
        }

        patch_smoother.update(make_ref(A));
        patch_smoother.apply(b, x);

        Vector r = b - A * x;
        Scalar norm_r = norm1(r);
        utopia_test_assert(norm_r < 1e-5);
    }
};

void patch_smoother() {
#ifdef UTOPIA_WITH_PETSC
    PatchSmootherTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(patch_smoother);
