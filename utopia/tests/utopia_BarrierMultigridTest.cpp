#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

#include "utopia_BarrierMultigrid.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_WITH_PETSC

namespace utopia {

    template <class Matrix, class Vector>
    class BarrierMultigridTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;
        using QuadraticFunction = utopia::QuadraticFunction<Matrix, Vector>;
        using Layout = typename Traits::Layout;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void run() {
            const static bool verbose = true;
            const static bool use_masks = false;
            int n_levels = 2;
            int n_coarse = 40;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);

            auto &&funs = ml_problem.get_functions();
            auto &&transfers = ml_problem.get_transfer();

            auto &&fun = funs.back();

            Vector x;
            fun->get_eq_constrains_values(x);

            Vector lower_bound(layout(x), -0.8), upper_bound(layout(x), 0.1);
            BoxConstraints<Vector> box(nullptr, make_ref(upper_bound));

            auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector>>();
            auto preconditioner = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            linear_solver->set_preconditioner(preconditioner);

            InputParameters params;
            params.set("use_coarse_space", false);
            params.set("debug", true);
            // params.set("show", true);
            params.set("barrier_parameter", 1);
            params.set("barrier_parameter_shrinking_factor", 0.3);
            // params.set("barrier_parameter_shrinking_factor", 0.9);
            params.set("min_barrier_parameter", 1e-10);
            // params.set("min_barrier_parameter", 1e-15);
            params.set("max_it", 20000);
            params.set("barrier_function_type", "BoundedLogBarrier");
            params.set("use_non_linear_residual", false);
            params.set("pre_smoothing_steps", 10);

            BarrierMultigrid<Matrix, Vector> mg(linear_solver);
            mg.verbose(verbose);
            mg.read(params);

            mg.set_box_constraints(box);
            mg.set_transfer_operators(transfers);
            mg.solve(*fun, x);

            rename("x", x);
            write("X.m", x);
        }
    };

    static void barrier_mg() {
#ifdef UTOPIA_WITH_PETSC
        BarrierMultigridTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC

        // #ifdef UTOPIA_WITH_TRILINOS
        //         BarrierMultigridTest<TpetraMatrixd, TpetraVectord>().run();
        // #endif  // UTOPIA_WITH_TRILINOS

        // #ifdef UTOPIA_WITH_BLAS
        //         BarrierMultigridTest<BlasMatrixd, BlasVectord>()
        //             .run();  // TODO(zulianp): : because blas is missing min operation ....
        // #endif               // UTOPIA_WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(barrier_mg);
}  // namespace utopia