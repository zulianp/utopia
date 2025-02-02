#include "test_problems/utopia_TestProblems.hpp"
#include "utopia.hpp"
#include "utopia_Rename.hpp"
#include "utopia_Testing.hpp"
#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {
    template <class Matrix>
    void assemble_laplacian_1D(const utopia::SizeType n, Matrix &m) {
        // n x n matrix with maximum 3 entries x row
        {
            Write<Matrix> w(m);
            Range r = row_range(m);

            // You can use set instead of add. [Warning] Petsc does not allow to mix add and set.
            for (SizeType i = r.begin(); i != r.end(); ++i) {
                if (i > 0) {
                    m.add(i, i - 1, -1.0);
                }

                if (i < n - 1) {
                    m.add(i, i + 1, -1.0);
                }

                m.add(i, i, 2.0);
            }
        }
    }

#ifdef UTOPIA_ENABLE_PETSC
    class PetscNonlinearSolverTest {
    public:
        using Traits = utopia::Traits<PetscVector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        Comm comm_;

        void run() {
            UTOPIA_RUN_TEST(petsc_ngs_test);
            UTOPIA_RUN_TEST(petsc_gss_newton_test);
            UTOPIA_RUN_TEST(petsc_newton_test);
            UTOPIA_RUN_TEST(petsc_sparse_newton_test_inexact);
            UTOPIA_RUN_TEST(petsc_newton_rosenbrock_test);
            UTOPIA_RUN_TEST(petsc_sparse_semismooth_newton_test);  // petsc 3.11.3 ERROR here
            UTOPIA_RUN_TEST(petsc_sparse_nonlinear_semismooth_newton_test);
            UTOPIA_RUN_TEST(petsc_direct_solver_newton_test);
            UTOPIA_RUN_TEST(petsc_newton_test_out_info);
            UTOPIA_RUN_TEST(petsc_sparse_newton_test);
            UTOPIA_RUN_TEST(petsc_newton_petsc_cg_test);
            UTOPIA_RUN_TEST(petsc_tr_rr_test);
            UTOPIA_RUN_TEST(petsc_snes_test);                // petsc 3.11.3 ERROR here
            UTOPIA_RUN_TEST(petsc_sparse_newton_snes_test);  // petsc 3.11.3 ERROR here
        }

        void petsc_ngs_test() {
            const int n = 50;
            PetscMatrix m;
            m.sparse(layout(comm_, Traits::decide(), Traits::decide(), n, n), 3, 2);
            assemble_laplacian_1D(n, m);
            // const Scalar ub = 100.0;
            const Scalar ub = 1.;
            const bool use_line_search = true;
            // const bool use_line_search = false;
            const int max_it = n * 60;
            const bool verbose = false;
            const int n_local_sweeps = 3;

            {
                Range r = row_range(m);
                Write<PetscMatrix> w(m);
                if (r.begin() == 0) {
                    m.set(0, 0, 1.);
                    m.set(0, 1, 0);
                }

                if (r.end() == n) {
                    m.set(n - 1, n - 1, 1.);
                    m.set(n - 1, n - 2, 0);
                }
            }

            PetscVector rhs(row_layout(m), 1.);
            rhs *= 1. / (n - 1);

            {
                // Creating test vector (alternative way see [assemble vector alternative], which might be easier for
                // beginners)
                Range r = range(rhs);
                Write<PetscVector> w(rhs);

                if (r.begin() == 0) {
                    rhs.set(0, 0);
                }

                if (r.end() == n) {
                    rhs.set(n - 1, 0.);
                }
            }

            PetscVector upper_bound(layout(rhs), ub);
            PetscVector solution(layout(rhs), 0.0);

            ProjectedGaussSeidel<PetscMatrix, PetscVector, -1> pgs;
            pgs.max_it(max_it);
            pgs.verbose(verbose);
            pgs.use_line_search(use_line_search);
            pgs.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

            Chrono c;
            c.start();

            pgs.solve(m, rhs, solution);

            c.stop();
            // if(mpi_world_rank() == 0) utopia::out() <<c << std::endl;

            PetscVector solution_u(layout(rhs), 0.0);
            ProjectedGaussSeidel<PetscMatrix, PetscVector, -1> pgs_u;
            pgs_u.verbose(verbose);
            pgs_u.n_local_sweeps(n_local_sweeps);
            pgs_u.use_line_search(use_line_search);
            pgs_u.use_symmetric_sweep(true);
            pgs_u.max_it(max_it);

            pgs_u.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

            c.start();

            pgs_u.solve(m, rhs, solution_u);

            c.stop();
            // if(mpi_world_rank() == 0) utopia::out() <<c << std::endl;

            Scalar diff = norm2(solution_u - solution);
            // Scalar res_norm = norm2(m * solution_u - rhs);

            // disp(res_norm);

            if (diff > 1e-5) {
                std::cerr << "[Error] different implementations of pgs gives different results, diff: " << diff
                          << std::endl;
            }

            utopia_test_assert(approxeq(solution_u, solution, 1e-5));

            // standard gs with MatSOR
            // GaussSeidel<PetscMatrix, PetscVector> gs;
            // gs.verbose(verbose);
            // gs.max_it(max_it);
            // // gs.sweeps(1);
            // solution.set(0.);

            // c.start();

            // gs.solve(m, rhs, solution);

            // c.stop();

            // if(mpi_world_rank() == 0) utopia::out() <<c << std::endl;

            // Scalar res_norm_ref = norm2(m * solution - rhs);

            // // disp(res_norm_ref);

            // if(diff > 1e-5) {
            // 	std::cerr << "[Error] different implementations of pgs gives different results, diff: " << diff <<
            // std::endl;
            // }

            // utopia_test_assert(approxeq(solution_u, solution, 1e-5));
        }

        void petsc_gss_newton_test() {
            using F = std::function<void(
                const PetscMatrix &, const PetscVector &, const PetscVector &, PetscVector &, PetscVector &)>;

            const int n = comm_.size() * 4;

            auto vl = layout(comm_, 4, n);
            PetscVector sol(vl, 0.0);
            PetscVector upbo(vl, 1.);
            PetscMatrix A;
            A.identity(square_matrix_layout(vl));
            PetscVector rhs(vl, 3.);

            PetscVector lambda, d;
            F f = [&lambda, &d, &upbo](const PetscMatrix &H,
                                       const PetscVector & /*g*/,
                                       const PetscVector &x,
                                       PetscVector &active,
                                       PetscVector &value) {
                lambda = (upbo - H * x);
                d = lambda + (x - upbo);

                Read<PetscVector> r_d(d);
                Read<PetscVector> r_u(upbo);
                Write<PetscVector> w_d(active);
                Write<PetscVector> w_v(value);

                auto rr = range(x);
                for (SizeType i = rr.begin(); i != rr.end(); i++) {
                    if (d.get(i) >= -1e-16) {
                        active.set(i, 1.0);
                        value.set(i, upbo.get(i));
                    } else {
                        active.set(i, 0.0);
                        value.set(i, 0.);
                    }
                }
            };

            auto linear_solver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();
            GenericSemismoothNewton<PetscMatrix, PetscVector, F> solver(f, linear_solver);

            solver.solve(A, rhs, sol);
        }

        void petsc_tr_rr_test() {
            // rosenbrock test
            if (comm_.size() == 1) {
                PetscVector x(serial_layout(10), 2);
                TestFunctionND_1<PetscMatrix, PetscVector> fun2(x.comm(), x.size());
                PetscVector expected(layout(x), 0.468919);

                InputParameters in;
                in.set("atol", 1e-10);
                in.set("rtol", 1e-10);
                in.set("stol", 1e-10);

                auto subproblem = std::make_shared<utopia::SteihaugToint<PetscMatrix, PetscVector>>();
                TrustRegion<PetscMatrix, PetscVector> tr_solver(subproblem);
                tr_solver.read(in);
                tr_solver.solve(fun2, x);

                x.set(2);
                trust_region_solve(fun2, x, Solver::steihaug_toint(), in);

                x.set(2);
                trust_region_solve(fun2, x, Solver::nash(), in);

                x.set(2);
                trust_region_solve(fun2, x, Solver::lanczos(), in);

                x.set(2);
                trust_region_solve(fun2, x, Solver::cgne(), in);
            }
        }

        void petsc_newton_test_out_info() {
            if (comm_.size() > 10) {
                return;
            }

            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);

            nlsolver.verbose(false);

            PetscVector x(layout(comm_, Traits::decide(), 10), 2);
            TestFunctionND_1<PetscMatrix, PetscVector> fun2(x.comm(), x.size());

            PetscVector expected(layout(x), 0.468919);
            nlsolver.solve(fun2, x);
            utopia_test_assert(approxeq(expected, x));
        }

        void petsc_sparse_newton_test() {
            if (comm_.size() > 10) {
                return;
            }

            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);
            nlsolver.verbose(false);

            const SizeType n = 10;
            SimpleQuadraticFunction<PetscMatrix, PetscVector> fun(n);

            PetscVector x(layout(comm_, Traits::decide(), n), 2.);
            PetscVector expected(layout(x), 0.0);

            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));
        }

        void petsc_sparse_newton_test_inexact() {
            if (comm_.size() > 10) {
                return;
            }

            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.verbose(false);
            nlsolver.forcing_strategy(InexactNewtonForcingStartegies::CAI);

            const SizeType n = 10;
            SimpleQuadraticFunction<PetscMatrix, PetscVector> fun(n);

            PetscVector x(layout(comm_, Traits::decide(), n), 2.);
            PetscVector expected(layout(x), 0.0);

            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));
        }

        void petsc_newton_test() {
            if (comm_.size() > 10) {
                return;
            }

            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);
            nlsolver.verbose(false);

            const SizeType n = 10;
            SimpleQuadraticFunction<PetscMatrix, PetscVector> fun(n);

            PetscVector x(layout(comm_, Traits::decide(), n), 2.);
            PetscVector expected(layout(x), 0.0);

            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));

            x.set(2.0);
            TestFunctionND_1<PetscMatrix, PetscVector> fun2(x.comm(), x.size());

            expected.set(0.468919);
            nlsolver.solve(fun2, x);
            utopia_test_assert(approxeq(expected, x));
        }

        void petsc_newton_rosenbrock_test() {
            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);
            nlsolver.rtol(1e-15);
            nlsolver.stol(1e-15);
            nlsolver.atol(1e-15);
            nlsolver.verbose(false);

            PetscVector expected_rosenbrock;
            PetscVector x0;

            if (comm_.size() <= 2) {
                expected_rosenbrock.values(serial_layout(2), 1.0);
                ExtendedRosenbrock21<PetscMatrix, PetscVector> r_generic_2d(comm_, local_size(expected_rosenbrock));
                x0.values(serial_layout(2), 2.0);
                nlsolver.solve(r_generic_2d, x0);
                utopia_test_assert(approxeq(expected_rosenbrock, x0));
            }

            if (comm_.size() <= 3) {
                expected_rosenbrock.values(serial_layout(3), 1.0);
                ExtendedRosenbrock21<PetscMatrix, PetscVector> r_generic_3d(comm_, local_size(expected_rosenbrock));
                x0.values(serial_layout(3), -2.0);
                nlsolver.solve(r_generic_3d, x0);
                utopia_test_assert(approxeq(expected_rosenbrock, x0));
            }

            if (comm_.size() <= 6) {
                expected_rosenbrock.values(serial_layout(6), 1.0);
                ExtendedRosenbrock21<PetscMatrix, PetscVector> r_generic_6d(comm_, local_size(expected_rosenbrock));
                x0.values(serial_layout(6), 2.0);
                nlsolver.solve(r_generic_6d, x0);
                utopia_test_assert(approxeq(expected_rosenbrock, x0));
            }
        }

        void petsc_sparse_semismooth_newton_test() {
            auto lsolver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();

            PetscMatrix A;
            PetscVector b, ub;

            bool verbose = false;

            SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> petsc_ss_newton(lsolver);
            SemismoothNewton<PetscMatrix, PetscVector, HOMEMADE> homemade_ss_newton(lsolver);

            petsc_ss_newton.verbose(verbose);
            homemade_ss_newton.verbose(verbose);

            InputParameters hm_params;
            // hm_params.set("use-adaptive-tol", true);
            homemade_ss_newton.read(hm_params);

            Poisson1D<PetscMatrix, PetscVector> ex2(_n, 2);
            PetscVector x_0 = ex2.initial_guess();
            ex2.hessian(x_0, A);
            ex2.get_rhs(b);
            ub = ex2.upper_bound();

            const Scalar scale_factor = 1;
            A *= scale_factor;
            b *= scale_factor;
            ub *= scale_factor;

            auto box = make_upper_bound_constraints(make_ref(ub));
            petsc_ss_newton.set_box_constraints(box);
            petsc_ss_newton.atol(1e-18);
            petsc_ss_newton.rtol(1e-15);
            petsc_ss_newton.stol(1e-16);
            petsc_ss_newton.max_it(400);
            petsc_ss_newton.solve(A, b, x_0);

            PetscVector hm_x_0(layout(b), 0.0);
            homemade_ss_newton.set_box_constraints(box);
            homemade_ss_newton.stol(1e-16);
            // homemade_ss_newton.verbose(true);
            homemade_ss_newton.solve(A, b, hm_x_0);

            // rename("x", x_0);
            // rename("y", hm_x_0);

            x_0 *= 1. / scale_factor;
            hm_x_0 *= 1. / scale_factor;

            // write("x_p.m", x_0);
            // write("x_u.m", hm_x_0);

            if (!approxeq(x_0, hm_x_0, 1e-14)) {
                PetscVector diff = hm_x_0 - x_0;

                Scalar norm_diff = norm2(diff) / Scalar(norm2(hm_x_0));
                utopia::out() << "norm_diff: " << std::to_string(norm_diff) << std::endl;
            }

            utopia_test_assert(approxeq(x_0, hm_x_0, 1e-8));
        }

        void petsc_sparse_nonlinear_semismooth_newton_test() {
            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();

            NonlinSemismoothNewton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.enable_differentiation_control(false);
            nlsolver.verbose(false);

            PetscVector upbo;

            Poisson1D<PetscMatrix, PetscVector> fun(_n, 3);
            PetscVector x_0 = fun.initial_guess();
            upbo = fun.upper_bound();

            auto box = make_upper_bound_constraints(make_ref(upbo));
            nlsolver.set_box_constraints(make_ref(box));

            nlsolver.solve(fun, x_0);
        }

        void petsc_direct_solver_newton_test() {
            auto lsolver = std::make_shared<Factorization<PetscMatrix, PetscVector>>();

#ifdef PETSC_HAVE_MUMPS
            lsolver->set_type(Solver::mumps(), Solver::lu_decomposition());
#endif  // PETSC_HAVE_MUMPS

            Newton<PetscMatrix, PetscVector> nlsolver(lsolver);
            nlsolver.verbose(false);

            SimpleQuadraticFunction<PetscMatrix, PetscVector> fun(_n);

            PetscVector x(layout(comm_, Traits::decide(), _n), 2.);
            PetscVector expected(layout(x), 0.0);

            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));

            auto lCG = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();
            nlsolver.set_linear_solver(lCG);
            x.set(2.);
            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));
        }

        void petsc_newton_petsc_cg_test() {
            using namespace std;

            // CG with diagonal preconditioner
            auto linear_solver = make_shared<KSPSolver<PetscMatrix, PetscVector>>();
            auto preconditioner = make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
            linear_solver->set_preconditioner(preconditioner);

            // Newton solver with cg linear solver
            Newton<PetscMatrix, PetscVector> newton_solver(linear_solver);
            newton_solver.verbose(false);

            const int n = 10;
            PetscVector actual(layout(comm_, Traits::decide(), n), 2.);
            PetscVector expected(layout(actual), 0.468919);

            TestFunctionND_1<PetscMatrix, PetscVector> fun(actual.comm(), n);

            newton_solver.solve(fun, actual);
            utopia_test_assert(approxeq(expected, actual));
        }

        void petsc_snes_test() {
            using namespace utopia;
            using namespace std;

            const static bool verbose = false;

            if (comm_.size() >= 10) {
                return;
            }

            auto linear_solver = make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();

            SNESSolver<PetscMatrix, PetscVector, PETSC> nonlinear_solver(linear_solver);
            nonlinear_solver.verbose(verbose);

            if (comm_.size() == 1) {
                Rosenbrock01<PetscMatrix, PetscVector> rosenbrock;
                PetscVector expected_rosenbrock(serial_layout(2), 1.0);
                PetscVector x0_ros(serial_layout(2), 1.5);

                nonlinear_solver.solve(rosenbrock, x0_ros);

                expected_rosenbrock -= x0_ros;
                Scalar diff_rb = norm2(expected_rosenbrock);
                utopia_test_assert(approxeq(diff_rb, 0., 1e-6));

                // std::cout<<"--------------------------------------------------- \n";
                auto cg_home = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();
                cg_home->verbose(verbose);

                SNESSolver<PetscMatrix, PetscVector, PETSC> nonlinear_solver2(cg_home);
                nonlinear_solver2.verbose(verbose);

                // reset IG
                x0_ros.values(serial_layout(2), 1.5);
                expected_rosenbrock.values(serial_layout(2), 1.0);
                nonlinear_solver2.solve(rosenbrock, x0_ros);

                expected_rosenbrock -= x0_ros;
                diff_rb = norm2(expected_rosenbrock);
                utopia_test_assert(approxeq(diff_rb, 0., 1e-6));

                // std::cout<<"------------------ utopia-precond test --------------------------------- \n";

                auto preconditioner = make_shared<InvDiagPreconditioner<PetscMatrix, PetscVector>>();
                cg_home->set_preconditioner(preconditioner);

                SNESSolver<PetscMatrix, PetscVector, PETSC> nonlinear_solver3(cg_home);
                nonlinear_solver3.verbose(verbose);
            }
        }

        void petsc_sparse_newton_snes_test() {
            if (comm_.size() > 1) {
                return;
            }

            auto lsolver = std::make_shared<BiCGStab<PetscMatrix, PetscVector>>();
            Newton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> nlsolver(lsolver);
            nlsolver.verbose(false);

            const SizeType n = 10;
            SimpleQuadraticFunction<PetscMatrix, PetscVector> fun(n);

            PetscVector x(layout(comm_, Traits::decide(), n), 2.);
            PetscVector expected(layout(x), 0.0);

            nlsolver.solve(fun, x);
            utopia_test_assert(approxeq(expected, x));

            if (comm_.size() == 1) {
                Rosenbrock01<PetscMatrix, PetscVector> rosenbrock;
                PetscVector expected_rosenbrock(serial_layout(2), 1.0);
                PetscVector x0_ros(serial_layout(2), 1.5);

                nlsolver.line_search_type("cp");
                nlsolver.line_search_order(3);
                nlsolver.max_it(1000);
                nlsolver.solve(rosenbrock, x0_ros);

                expected_rosenbrock -= x0_ros;
                Scalar diff_rb = norm2(expected_rosenbrock);
                UTOPIA_UNUSED(diff_rb);
                utopia_test_assert(approxeq(diff_rb, 0., 1e-6));
            }
        }

        PetscNonlinearSolverTest() = default;

    private:
        int _n{100};
    };

#endif  // UTOPIA_ENABLE_PETSC

    static void petsc_nonlinear() {
#ifdef UTOPIA_ENABLE_PETSC
        PetscNonlinearSolverTest().run();
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(petsc_nonlinear);
}  // namespace utopia
