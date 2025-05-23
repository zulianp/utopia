#include "utopia_Testing.hpp"

#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_PETSC

#include "utopia.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

#include "utopia_BarrierMultigrid.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_LogBarrierQPMultigrid.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_ENABLE_PETSC

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

        static std::unique_ptr<BarrierMultigrid<Matrix, Vector>> create_barrier_mg(const int n_levels,
                                                                                   const bool algebraic,
                                                                                   bool verbose) {
            auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
            auto preconditioner = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            linear_solver->set_preconditioner(preconditioner);
            linear_solver->max_it(10000);

            // auto linear_solver = std::make_shared<Factorization<Matrix, Vector>>();

            InputParameters params;
            params.set("use_coarse_space", true);
            // params.set("debug", true);
            params.set("verbose", verbose);
            params.set("barrier_parameter", 1);
            params.set("barrier_parameter_shrinking_factor", 0.1);
            params.set("min_barrier_parameter", 1e-10);
            params.set("max_it", 40);
            params.set("barrier_function_type", "BoundedLogBarrier");
            // params.set("barrier_function_type", "LogBarrier");

            params.set("use_non_linear_residual", false);
            params.set("pre_smoothing_steps", 5);
            params.set("post_smoothing_steps", 5);
            params.set("atol", 1e-7);
            params.set("rtol", 1e-8);
            params.set("mg_steps", 1);
            params.set("keep_initial_coarse_spaces", true);
            params.set("amg_n_coarse_spaces", n_levels - 1);
            params.set("non_smooth_projection_weight", 0.9);

            auto mg = utopia::make_unique<BarrierMultigrid<Matrix, Vector>>(linear_solver);
            mg->verbose(verbose);
            mg->read(params);

            // Use external linear smoothing (otherwise internally uses Jacobi)
            // mg->set_linear_smoother(std::make_shared<ILU<Matrix, Vector>>());

#ifdef UTOPIA_ENABLE_PETSC
            if (algebraic) {
                auto agg = std::make_shared<Agglomerate<Matrix>>();

                InputParameters agg_params;
                agg_params.set("bmax", 4);
                agg->read(agg_params);

                agg->verbose(false);
                mg->set_agglomerator(agg);
            }
#endif

            return mg;
        }

        void test_ml_problem() {
            const static bool verbose = false;
            const static bool use_masks = false;
            int n_levels = 5;
            int n_coarse = 501;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);

            auto &&funs = ml_problem.get_functions();
            auto &&transfers = ml_problem.get_transfer();

            auto &&fun = funs.back();

            Vector x;
            fun->get_eq_constrains_values(x);

            Vector lower_bound(layout(x), -0.8), upper_bound(layout(x), 0.8);
            BoxConstraints<Vector> box(nullptr, make_ref(upper_bound));

            // bool algebraic = Traits::Backend == PETSC;
            bool algebraic = false;

            auto mg = create_barrier_mg(n_levels, algebraic, verbose);
            mg->set_box_constraints(box);

            if (!algebraic) {
                mg->set_transfer_operators(transfers);
            }

            mg->solve(*fun, x);

            // if (Traits::Backend == PETSC) {
            //     rename("x", x);
            //     write("X.m", x);
            // }
        }

#ifdef UTOPIA_ENABLE_PETSC
        void test_qp_problem() {
            //////////////////////////////////////////////////////////////////////////////////////////
            // Problem set-up
            const static bool verbose = false;
            const static bool use_masks = false;
            int n_levels = 5;
            int n_coarse = 501;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);

            Matrix A;
            Vector x, b;

            auto &&funs = ml_problem.get_functions();
            auto &&fun = funs.back();
            fun->get_eq_constrains_values(x);
            fun->gradient(x, b);
            fun->hessian(x, A);

            // Negative gradient!
            b = -b;

            Vector lower_bound(layout(x), -0.8), upper_bound(layout(x), 0.8);
            BoxConstraints<Vector> box(nullptr, make_ref(upper_bound));

            //////////////////////////////////////////////////////////////////////////////////////////
            // Solver set-up
            auto linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
            auto preconditioner = std::make_shared<InvDiagPreconditioner<Matrix, Vector>>();
            linear_solver->set_preconditioner(preconditioner);
            linear_solver->max_it(10000);

            LogBarrierQPMultigrid<Matrix, Vector> mg(linear_solver);

            mg.set_box_constraints(box);
            mg.verbose(verbose);

            // Solves
            mg.solve(A, b, x);

            if (Traits::Backend == PETSC) {
                rename("x", x);
                write("X_qp.m", x);
            }
        }
#endif  // UTOPIA_ENABLE_PETSC

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(test_ml_problem);
#ifdef UTOPIA_ENABLE_PETSC
            UTOPIA_RUN_TEST(test_qp_problem);
#endif
        }
    };

    static void barrier_mg() {
#ifdef UTOPIA_ENABLE_PETSC
        BarrierMultigridTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_ENABLE_PETSC

        // #ifdef UTOPIA_ENABLE_TRILINOS
        //         BarrierMultigridTest<TpetraMatrixd, TpetraVectord>().run();
        // #endif  // UTOPIA_ENABLE_TRILINOS

        // #ifdef UTOPIA_ENABLE_BLAS
        //         BarrierMultigridTest<BlasMatrixd, BlasVectord>()
        //             .run();  // TODO(zulianp): : because blas is missing min operation ....
        // #endif               // UTOPIA_ENABLE_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(barrier_mg);
}  // namespace utopia

#endif
