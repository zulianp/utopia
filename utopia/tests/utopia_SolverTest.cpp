#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_TestProblems.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_MSSolver.hpp"
#include "utopia_SPStaticCondensationKrylov.hpp"

namespace utopia
{

    template<class Matrix, class Vector, class Scalar>
    class SolverTest {
    public:
        static void print_backend_info()
        {
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void run()
        {
            print_backend_info();
            UTOPIA_RUN_TEST(ms_solver);
            UTOPIA_RUN_TEST(newton_cg_test);
            UTOPIA_RUN_TEST(grad_descent_test);
            UTOPIA_RUN_TEST(tr_test);
            UTOPIA_RUN_TEST(ls_test);
            UTOPIA_RUN_TEST(nl_solve_test);
            UTOPIA_RUN_TEST(dogleg_test);
            UTOPIA_RUN_TEST(st_cg_test);
            UTOPIA_RUN_TEST(precond_st_cg_test);
        }

        class EmptyLSFun : public LeastSquaresFunction<Matrix, Vector> {
        public:
            bool residual(const Vector &, Vector &) const { return true; }
            bool jacobian(const Vector &, Matrix &) const { return true; }
            bool value(const Vector &, Scalar &) const { return true; }
            bool update(const Vector &) { return true; }
        };


        void ms_solver()
        {
            const int n = 20;
            // Rosenbrock01<Matrix, Vector> fun;
            TestFunctionND_1<Matrix, Vector> fun(n);


            // SimpleQuadraticFunction<Matrix, Vector> fun;
            // Rastrigin<Matrix, Vector> fun;
            // Woods14<Matrix, Vector> fun;

            Vector x = values(n, 2.0);

            MSSolver<Matrix, Vector> solver(std::make_shared<ConjugateGradient<Matrix, Vector>>());
            solver.set_norm_type(MSSolver<Matrix, Vector>::A_SQUARED_NORM);
            // solver.set_norm_type(MSSolver<Matrix, Vector>::A_NORM);
            // solver.verbose(true);
            solver.solve(fun, x);
            // disp(x);
        }

        void ls_normal_eq()
        {
            LeastSquaresNewton<Matrix, Vector> newton(std::make_shared<ConjugateGradient<Matrix, Vector>>());
            auto ls_strat  = std::make_shared<utopia::Backtracking<Vector> >();
            newton.set_line_search_strategy(ls_strat);

            EmptyLSFun fun;
            Vector x0;
            newton.solve(fun, x0);
        }

        void st_cg_test()
        {
            SteihaugToint<Matrix, Vector, HOMEMADE> cg;
            cg.rtol(1e-7);
            cg.atol(1e-7);
            cg.max_it(_n);
            cg.verbose(false);

            Matrix A = sparse(_n, _n, 3);
            assemble_symmetric_laplacian_1D(A, true);

            Vector rhs = values(_n, 975.9);

            {
                auto r = range(rhs);
                Write<Vector> w(rhs);
                if(r.inside(0)) rhs.set(0, 0.0);
                if(r.inside(_n - 1)) rhs.set(_n-1, 0.0);
            }

            Vector x = zeros(size(rhs));

            cg.solve(A, rhs, x);
            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }

        void precond_st_cg_test()
        {
            SteihaugToint<Matrix, Vector, HOMEMADE> cg;
            cg.rtol(1e-7);
            cg.atol(1e-7);
            cg.max_it(_n);
            cg.verbose(false);
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
            // cg.set_preconditioner(std::make_shared<IdentityPreconditioner<Matrix, Vector> >());


            Matrix A = sparse(_n, _n, 3);
            assemble_symmetric_laplacian_1D(A, true);

            Vector rhs = values(_n, 975.9);

            {
                auto r = range(rhs);
                Write<Vector> w(rhs);
                if(r.inside(0)) rhs.set(0, 0.0);
                if(r.inside(_n-1)) rhs.set(_n-1, 0.0);
            }

            Vector x = zeros(size(rhs));

            cg.solve(A, rhs, x);
            utopia_test_assert(approxeq(rhs, A * x, 1e-5));
        }

        void nl_solve_test()
        {
            //! [NL solve example]

            //set-up problem
            int n = 10;
            Vector actual = values(n, 2.);
            TestFunctionND_1<Matrix, Vector> fun(n);

            InputParameters in;
            in.set("atol", 1e-11);
            in.set("verbose", false);


            //solve problem
            line_search_solve(fun, actual, Solver::backtracking(), in);

            //test outcome...
            Vector expected = values(n, 0.468919);
            utopia_test_assert(approxeq(expected, actual));
            //! [NL solve example]
        }

        void grad_descent_test()
        {
            //set-up problem
            int n = 10;
            Vector actual = values(n, 1.0);
            TestFunctionND_1<Matrix, Vector> fun(n);

            auto solver = GradientDescent<Vector>();
            solver.dumping_parameter(0.05);
            solver.solve(fun, actual);

            //test outcome...
            Vector expected = values(n, 0.468919);
            utopia_test_assert(approxeq(expected, actual));
        }


        void newton_cg_test()
        {
            //! [Newton CG example]
            using namespace std;

            //CG with diagonal preconditioner
            auto linear_solver  = make_shared< ConjugateGradient<Matrix, Vector> >();
            auto preconditioner = make_shared< InvDiagPreconditioner<Matrix, Vector> >();
            linear_solver->set_preconditioner(preconditioner);

            //Newton solver with cg linear solver
            Newton<Matrix, Vector> newton_solver(linear_solver);

            const int n = 10;
            Vector actual   = values(n, 2.);
            Vector expected = values(n, 0.468919);

            TestFunctionND_1<Matrix, Vector> fun(n);

            newton_solver.solve(fun, actual);

            //Check if the result is what we expected
            utopia_test_assert(approxeq(expected, actual));
            //! [Newton CG example]
        }

        void tr_test()
        {
            // rosenbrock test
            if(mpi_world_size() == 1)
            {
                Vector x = values(10, 2);
                TestFunctionND_1<Matrix, Vector> fun2(x.size());
                Vector expected = values(x.size(), 0.468919);

                Vector x_w1  = values(4, 10);
                Vector expected_woods = values(4, 1);
                Woods14<Matrix, Vector> fun_woods;
                {
                    Write<Vector> w1(x_w1);
                    x_w1.set(0, -3);
                    x_w1.set(1, -1);
                    x_w1.set(2, -3);
                    x_w1.set(3, -1);
                }

                InputParameters in;
                in.set("atol", 1e-10);
                in.set("rtol", 1e-10);
                in.set("stol", 1e-10);
                in.set("verbose", false);

                x = values(10, 2);
                trust_region_solve(fun2, x, Solver::steihaug_toint(), in);
                utopia_test_assert(approxeq(expected, x));

                x = values(10, 2);
                trust_region_solve(fun2, x, Solver::dogleg(), in);
                utopia_test_assert(approxeq(expected, x));

                x = values(10, 2);
                trust_region_solve(fun2, x, Solver::cauchypoint(), in);
                utopia_test_assert(approxeq(expected, x));

                Vector expected_rosenbrock = values(2, 1.0);
                Rosenbrock01<Matrix, Vector> rosenbrock;
                Vector x0 = values(2, 2.0);

                x0 = values(2, 2.0);
                in.set("atol", 1e-13);
                in.set("rtol", 1e-17);
                trust_region_solve(rosenbrock, x0, Solver::steihaug_toint(), in);

                auto diff_norm = norm_infty(expected_rosenbrock - x0);


                if(diff_norm > 1e-11) {
                    utopia_error("tr_test: Solver::steihaug_toint() with rosenbrock is failing");
                }
            }
        }

        void ls_test()
        {
            if(mpi_world_size() == 1) {
                Vector x1 = values(10, 2);
                Vector x2 = values(10, 2);
                TestFunctionND_1<Matrix, Vector> fun2(x1.size());

                Vector expected = values(x1.size(), 0.468919);

                InputParameters params;
                params.set("atol", 1e-11);
                params.set("rtol", 1e-11);
                params.set("stol", 1e-11);
                params.set("verbose", false);


                auto lsolver = std::make_shared< ConjugateGradient<Matrix, Vector, HOMEMADE> >();
                Newton<Matrix, Vector> nlsolver1(lsolver);
                Newton<Matrix, Vector> nlsolver2(lsolver);


                auto strategy_sbc = std::make_shared<utopia::SimpleBacktracking<Vector> >();
                auto strategy_bc  = std::make_shared<utopia::Backtracking<Vector, HOMEMADE> >();


                nlsolver1.set_line_search_strategy(strategy_sbc);
                nlsolver2.set_line_search_strategy(strategy_bc);

                nlsolver1.read(params);
                nlsolver2.read(params);


                // nlsolver1.solve(fun2, x1);
                // nlsolver2.solve(fun2, x2);

                

                // Woods function test
                Vector x_w1  = values(4, 10);
                Vector x_w2  = values(4, 10);
                Vector expected_woods = values(4, 1);
                {
                    Write<Vector> w1(x_w1);
                    Write<Vector> w2(x_w2);
                    x_w1.set(0, -3);    x_w2.set(0, -3);
                    x_w1.set(1, -1);    x_w2.set(1, -1);
                    x_w1.set(2, -3);    x_w2.set(2, -3);
                    x_w1.set(3, -1);    x_w2.set(3, -1);
                }

                Woods14<Matrix, Vector> fun_woods;
                nlsolver1.solve(fun_woods, x_w1);
                nlsolver2.solve(fun_woods, x_w2);


                utopia_test_assert(approxeq(expected_woods, x_w1));
                utopia_test_assert(approxeq(expected_woods, x_w2));

                // rastrigin function test - convergence to local minimum
                Rastrigin<Matrix, Vector> fun_rastrigin(2);
                Vector x_r1 = values(2, 1), x_r2 = values(2, 1), expected_rastrigin = values(2, 1);
                {
                    Write<Vector> w1(x_r1);
                    Write<Vector> w2(x_r2);
                    Write<Vector> w3(expected_rastrigin);
                    x_r1.set(0, -5.12); x_r2.set(0, -5.12);  expected_rastrigin.set(0, -4.97469);
                    x_r1.set(1, 5.12); x_r2.set(1, 5.12);  expected_rastrigin.set(1, 4.97469);
                }

                nlsolver1.solve(fun_rastrigin, x_r2);
                nlsolver2.solve(fun_rastrigin, x_r2);

                // rosenbrock test
                Vector expected_rosenbrock = values(2, 1);
                Rosenbrock01<Matrix, Vector> rosenbrock_fun;

                Vector x01 = values(2, 2.0), x02 = values(2, 2.0);
                nlsolver1.solve(rosenbrock_fun, x01);
                nlsolver2.solve(rosenbrock_fun, x02);


                utopia_test_assert(approxeq(expected_rosenbrock, x01));
                utopia_test_assert(approxeq(expected_rosenbrock, x02));
            }
        }


        void dogleg_test()
        {
            // rosenbrock test
            if(mpi_world_size() == 1)
            {
                Rosenbrock01<Matrix, Vector> rosenbrock;
                Vector expected_rosenbrock = values(2, 1);

                auto cg = std::make_shared<ConjugateGradient<Matrix, Vector> >();
                auto dogleg = std::make_shared<Dogleg<Matrix, Vector> >(cg);

                Vector x0 = values(2, 2.0);

                TrustRegion<Matrix, Vector> tr_solver(dogleg);
                tr_solver.verbose(false);
                tr_solver.max_it(100);
                tr_solver.solve(rosenbrock, x0);

                utopia_test_assert(approxeq(expected_rosenbrock, x0));
            }
        }


        SolverTest()
        : _n(10) { }

    private:
        int _n;

    };

    template<class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
    class MSSolverTest {
    public:

        void run()
        {
            //UTOPIA_UNIT_TEST_BEGIN("MSSolverTest");
            // UTOPIA_RUN_TEST(convex_hull_2);
            UTOPIA_RUN_TEST(convex_hull_4);
            // UTOPIA_RUN_TEST(convex_hull_8);
            //UTOPIA_UNIT_TEST_END("MSSolverTest");
        }

        void convex_hull(const int convex_hull_n_gradients)
        {
            const bool verbose = true;

            if(verbose) { std::cout << "Rastrigin:" << std::endl; }

            Rastrigin<GlobalMatrix, GlobalVector> fun1;
            aux_convex_hull(20, fun1, convex_hull_n_gradients);

            if(mpi_world_size() == 1) {
                //FIXME seems to fail for this function
                // if(verbose) { std::cout << "Rosenbrock01:" << std::endl; }

                // Rosenbrock01<GlobalMatrix, GlobalVector> fun2;
                // aux_convex_hull(2, fun2, convex_hull_n_gradients);

                if(verbose) { std::cout << "Woods14:" << std::endl; }

                Woods14<GlobalMatrix, GlobalVector> fun3;
                aux_convex_hull(4, fun3, convex_hull_n_gradients);
            }

            int n = 20;
            TestFunctionND_1<GlobalMatrix, GlobalVector> fun4(n);
            aux_convex_hull(n, fun4, convex_hull_n_gradients);
        }

        void convex_hull_2()
        {
            convex_hull(2);
        }

        void convex_hull_4()
        {
            convex_hull(4);
        }

        void convex_hull_8()
        {
            convex_hull(8);
        }

        void aux_convex_hull(
            const int n,
            Function<GlobalMatrix, GlobalVector> &fun,
            const int convex_hull_n_gradients)
        {
            using ConvexHullSolver = utopia::MSConvexHullSolver<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector>;

            GlobalVector x = values(n, 2.0);

            MSSolver<GlobalMatrix, GlobalVector> solver(std::make_shared<ConjugateGradient<GlobalMatrix, GlobalVector, HOMEMADE>>());
            // solver.set_norm_type(MSSolver<GlobalMatrix, GlobalVector>::A_SQUARED_NORM);
            // solver.set_norm_type(MSSolver<GlobalMatrix, GlobalVector>::A_NORM);

            solver.set_convex_hull_n_gradients(convex_hull_n_gradients);
            solver.set_convex_hull_solver(std::make_shared<ConvexHullSolver>());

            solver.verbose(true);
            // solver.atol(1e-10);
            solver.solve(fun, x);
        }

    };

    static void solvers()
    {

#ifdef WITH_BLAS
        // SolverTest<BlasMatrixd, BlasVectord, double>().run();
        //FIXME this fails for some reason
        // MSSolverTest<Matrixd, Vectord, Matrixd, Vectord>().run();
#endif //WITH_BLAS


#ifdef WITH_PETSC
        SolverTest<PetscMatrix, PetscVector, PetscScalar>().run();

#ifdef WITH_BLAS
        //FIXME this fails for some reason
        // MSSolverTest<PetscMatrix, PetscVector, Matrixd, Vectord>().run();
#endif //WITH_BLAS
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(solvers);
}
