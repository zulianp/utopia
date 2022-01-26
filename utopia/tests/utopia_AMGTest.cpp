#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_ILU.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_AlgebraicMultigrid.hpp"
#include "utopia_ILUDecompose.hpp"

#include "utopia_InputParameters.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_ILUDecompose.hpp"
#endif  // UTOPIA_WITH_PETSC

using namespace utopia;

template <class Matrix, class Vector>
class AMGTest {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using Comm = typename Traits::Communicator;

    SizeType n_levels{6};
    SizeType n_dofs = 100;

    void run() {
        // UTOPIA_RUN_TEST(test_amg_add_corr);
        UTOPIA_RUN_TEST(test_amg);
    }

    void test_amg() {
        // auto smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
        auto smoother = std::make_shared<ILU<Matrix, Vector>>();
        // auto coarse_solver = std::make_shared<ILU<Matrix, Vector>>();
        auto coarse_solver = std::make_shared<Factorization<Matrix, Vector>>();
        auto agglomerator = std::make_shared<Agglomerate<Matrix>>();

        AlgebraicMultigrid<Matrix, Vector> amg(smoother, coarse_solver, agglomerator);
        InputParameters params;
        params.set("n_levels", n_levels);
        amg.read(params);

        test_1D("amg_interp", amg);
    }

    void test_amg_add_corr() {
        // auto smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
        auto smoother = std::make_shared<ILU<Matrix, Vector>>();
        // auto coarse_solver = std::make_shared<ILU<Matrix, Vector>>();
        auto coarse_solver = std::make_shared<Factorization<Matrix, Vector>>();
        auto agglomerator = std::make_shared<Agglomerate<Matrix>>();
        agglomerator->use_additive_correction(true);

        AlgebraicMultigrid<Matrix, Vector> amg(smoother, coarse_solver, agglomerator);
        InputParameters params;
        params.set("n_levels", n_levels);
        amg.read(params);

        test_1D("add_corr", amg);
    }

    template <class AMG>
    void test_1D(const std::string &name, AMG &amg) {
        UTOPIA_UNUSED(name);

        auto comm = Comm::get_default();

        auto vl = layout(comm, n_dofs, n_dofs * comm.size());
        Matrix A;
        A.sparse(square_matrix_layout(vl), 3, 2);

        Vector x(vl, 0.0), b(vl, 2.0);

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

        // auto params = param_list(param(
        //     "measure_residual",
        //     param_list(param("type", "MeasureResidualComponents"), param("verbose", true), param("block_size", 2))));
        // amg.read(params);

        amg.verbose(false);
        amg.atol(1e-6);
        amg.rtol(1e-7);
        amg.stol(1e-7);
        amg.max_it(100);

        // amg.solve(A, b, x);

        amg.update(make_ref(A));

        // for (SizeType i = 0; i < n_levels; ++i) {
        //     auto &Ai = amg.algorithm().level(i).A();
        //     rename("a" + std::to_string(i) + "_" + name, const_cast<Matrix &>(Ai));
        //     write("A" + std::to_string(i) + "_" + name + ".m", Ai);
        // }

        amg.apply(b, x);

        Vector r = b - A * x;
        Scalar norm_r = norm1(r);
        utopia_test_assert(norm_r < 1e-5);
    }
};

void amg() {
#ifdef UTOPIA_WITH_PETSC
    AMGTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(amg);
