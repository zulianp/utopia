#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_PatchSmoother.hpp"

using namespace utopia;

template <class Matrix, class Vector>
class PatchSmootherTest {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using Comm = typename Traits::Communicator;

    SizeType n_levels{6};
    SizeType n_dofs = 100;

    void run() { UTOPIA_RUN_TEST(test_linear_patch_smoother); }

    void test_linear_patch_smoother() {
        {
            SOR<Matrix, Vector> gs;
            gs.verbose(true);
            gs.max_it(50000);

            test_1D("smoother_gs", gs);
        }

        {
            auto factorization = std::make_shared<Factorization<Matrix, Vector>>();

            PatchSmoother<Matrix> patch_smoother;
            patch_smoother.set_patch_solver(factorization);

            test_1D("patch_smoother", patch_smoother);
        }
    }

    template <class PatchSmoother>
    void test_1D(const std::string &name, PatchSmoother &patch_smoother) {
        UTOPIA_UNUSED(name);

        auto comm = Comm::get_default();

        auto vl = layout(comm, n_dofs, n_dofs * comm.size());
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
