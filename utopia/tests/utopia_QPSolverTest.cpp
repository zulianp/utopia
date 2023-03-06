#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_polymorphic_QPSolver.hpp"

#include "utopia_MonotoneAlgebraicMultigrid.hpp"
#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_MultigridQR.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_LogBarrierFunction.hpp"
#include "utopia_LogBarrierQPSolver.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_PrimalInteriorPointSolver_impl.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_BDDLinearSolver.hpp"
#include "utopia_petsc_BDDOperator.hpp"
#include "utopia_petsc_BDDQPSolver.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_ENABLE_PETSC

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolverTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        template <class QPSolver>
        void run_qp_solver(QPSolver &qp_solver) const {
            QPSolverTestProblem<Matrix, Vector>::run(n, verbose, qp_solver);
        }

        void pg_test() const {
            ProjectedGradient<Matrix, Vector> pg;
            run_qp_solver(pg);
        }

        void pcg_test() const {
            ProjectedConjugateGradient<Matrix, Vector> pcg;
            run_qp_solver(pcg);
        }

        void ngs_test() const {
            ProjectedGaussSeidel<Matrix, Vector> pgs;
            run_qp_solver(pgs);
        }

        void nblockgs_test() const {
            InputParameters params;
            params.set("block_size", 2);

            ProjectedGaussSeidel<Matrix, Vector> pgs;
            pgs.read(params);

            run_qp_solver(pgs);
        }

        static void create_symm_lapl_test_data(Comm &comm,
                                               Matrix &A,
                                               Vector &b,
                                               BoxConstraints<Vector> &box,
                                               SizeType n = 100,
                                               const bool boundary_conds = true) {
            A.sparse(layout(comm, Traits::decide(), Traits::decide(), n, n), 3, 2);
            assemble_symmetric_laplacian_1D(A, boundary_conds);

            auto h = 1. / (n - 1.);
            A = 1. / h * A;

            {
                Range r = row_range(A);
                const SizeType r_begin = r.begin();
                const SizeType r_end = r.end();

                Write<Matrix> w(A);
                if (r_begin == SizeType(0)) {
                    A.set(0, 0, 1.);
                }

                if (r_end == n) {
                    A.set(n - 1, n - 1, 1.);
                }
            }

            b.values(row_layout(A), 50.0);

            {
                Range row_range = range(b);
                const SizeType r_begin = row_range.begin();
                const SizeType r_end = row_range.end();

                Write<Vector> w(b);

                for (SizeType r = r_begin; r != r_end; ++r) {
                    if (r >= n / 2.) {
                        b.set(r, -50.0);
                    }
                    if (r == 0) {
                        b.set(r, 0);
                    }

                    if (r == (n - 1)) {
                        b.set(r, 0);
                    }
                }
            }

            b = h * b;

            auto lb = std::make_shared<Vector>(row_layout(A), -0.5);
            auto ub = std::make_shared<Vector>(row_layout(A), 0.5);

            box = make_box_constaints(lb, ub);
        }

        void log_barrier_test() {
            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box);

            // b *= 0.5;
            // box.lower_bound() = nullptr;
            // box.upper_bound() = nullptr;
            QuadraticFunction<Matrix, Vector> fun(make_ref(A), make_ref(b));

            InputParameters params;
            // params.set("verbose", true);
            params.set("barrier_parameter", 1e-5);
            params.set("barrier_parameter_shrinking_factor", 0.7);

            LogBarrierFunction<Matrix, Vector> barrier(make_ref(fun),
                                                       std::make_shared<LogBarrier<Matrix, Vector>>(make_ref(box)));
            barrier.read(params);

            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            // cg.verbose(true);
            cg.max_it(20);

            Newton<Matrix, Vector> newton(make_ref(cg));
            // newton.verbose(true);

            Vector x(layout(b));

            // Linear solve first to get closer to solution
            cg.solve(A, b, x);
            barrier.project_onto_feasibile_region(x);
            newton.solve(barrier, x);
        }

        void log_barrier_qp_solver_test(const std::string &barrier_function_type, const bool verbose) {
            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            SizeType n = 100;
            create_symm_lapl_test_data(comm, A, b, box, n);
            box.lower_bound() = nullptr;
            // box.upper_bound() = nullptr;

            InputParameters params;
            // params.set("verbose", true);
            params.set("barrier_parameter", 1e-2);
            params.set("barrier_thickness", 0.01);
            params.set("barrier_parameter_shrinking_factor", 0.1);
            params.set("min_barrier_parameter", 1e-5);
            params.set("verbose", verbose);
            params.set("function_type", barrier_function_type);
            params.set("max_it", 200);
            params.set("enable_line_search", true);

            LogBarrierQPSolver<Matrix, Vector> solver;

            ConjugateGradient<Matrix, Vector, HOMEMADE> linear_solver;
            linear_solver.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            linear_solver.max_it(100);
            solver.set_linear_solver(make_ref(linear_solver));

            solver.set_box_constraints(box);
            solver.read(params);

            Vector x(layout(b), 0.);
            utopia_test_assert(solver.solve(A, b, x));

            if (Traits::Backend == PETSC) {
                rename("x", x);
                write("X.m", x);

                if (box.has_lower_bound()) {
                    rename("lb", *box.lower_bound());
                    write("LB.m", *box.lower_bound());
                }

                if (box.has_upper_bound()) {
                    rename("ub", *box.upper_bound());
                    write("UB.m", *box.upper_bound());
                }

                rename("a", A);
                write("A.m", A);

                rename("b", b);
                write("B.m", b);
            }
        }

        void log_barrier_qp_solver_test() {
            // log_barrier_qp_solver_test("LogBarrierFunction", true);
            log_barrier_qp_solver_test("LogBarrierFunction", false);
            log_barrier_qp_solver_test("BoundedLogBarrierFunction", false);

            // Utopia::Abort("BYE");
        }

        void interior_point_qp_solver_test() {
            auto &&comm = Comm::get_default();

            SizeType n = 400;
            if (Traits::Backend == BLAS) {
                n = 30;
            }

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box, n);
            box.lower_bound() = nullptr;

            Vector selector(layout(b), 1.);

            if (true)  // Enable partial selection
            {
                Scalar start = 0.2;
                Scalar end = 0.25;
                auto selector_view = view_device(selector);
                parallel_for(
                    range_device(selector), UTOPIA_LAMBDA(const SizeType i) {
                        Scalar x = i * 1. / (n - 1);
                        bool val = x >= start && x <= end;
                        selector_view.set(i, val);
                    });
            }

            // bool verbose = Traits::Backend == PETSC;
            bool verbose = false;

            PrimalInteriorPointSolver<Matrix, Vector> solver;

            InputParameters params;
            params.set("min_val", 1e-18);
            params.set("debug", verbose);
            solver.read(params);

            ConjugateGradient<Matrix, Vector, HOMEMADE> linear_solver;
            linear_solver.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            linear_solver.max_it(10000);
            solver.set_linear_solver(make_ref(linear_solver));

            solver.set_box_constraints(box);
            solver.set_selection(make_ref(selector));
            solver.verbose(verbose);

            Vector x(layout(b), 0.);
            solver.solve(A, b, x);

            if (verbose) {
                rename("x", x);
                write("IP_X.m", x);
            }
        }

        void MPRGP_test() const {
            MPRGP<Matrix, Vector> qp_solver;
            run_qp_solver(qp_solver);

            auto &&comm = Comm::get_default();

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;
            create_symm_lapl_test_data(comm, A, b, box);

            Vector x = 0 * b;

            qp_solver.set_box_constraints(box);
            qp_solver.verbose(false);
            qp_solver.max_it(n * 2);
            qp_solver.set_eig_comp_tol(1e-1);
            qp_solver.solve(A, b, x);
        }

        void MG_QR_test() {
            bool verbose = false;

            Vector rhs, x;
            Vector upper_bound, lower_bound;
            Matrix A, R, Q, Ih_fine, Rot;
            Matrix Ih1, Ih0;

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/forQR/b", rhs);
            read(data_path + "/forQR/x", x);
            read(data_path + "/forQR/A", A);
            read(data_path + "/forQR/Q", Q);
            read(data_path + "/forQR/R", R);
            read(data_path + "/forQR/Rot", Rot);
            read(data_path + "/forQR/ub", upper_bound);
            read(data_path + "/forQR/lb", lower_bound);

            read(data_path + "/forQR/Ih", Ih_fine);
            read(data_path + "/forQR/I2h", Ih1);
            read(data_path + "/forQR/I3h", Ih0);

            x.set(1);

            auto num_levels = 3;

            // chop_abs(Q, 1e-7);

            // std::cout<<"A: "<< local_size(A).get(0) << "  \n";
            // std::cout<<"Q: "<< local_size(Q).get(0) << "  \n";
            // std::cout<<"R: "<< local_size(R).get(0) << ","<<local_size(R).get(1) <<
            // "  \n"; std::cout<<"I: "<< local_size(Ih_fine).get(0) << "  \n";
            // std::cout<<"rhs: "<< local_size(rhs).get(0) << "  \n";

            R = transpose(R);

            // version 1
            Matrix QtAQ = transpose(Q) * Rot * A * Rot * Q;
            Matrix QtIh = transpose(Q) * Rot * Ih_fine;
            Vector Qtrhs = transpose(Q) * Rot * rhs;
            Vector Qtx = transpose(Q) * Rot * x;

            // Matrix QtAQ  = Rot*A*Rot;
            // Matrix QtIh  = Rot* Ih_fine;
            // Vector Qtrhs = Rot *rhs;
            // Vector Qtx   = Rot *x;

            auto fine_smoother = std::make_shared<ProjectedGaussSeidelQR<Matrix, Vector>>();
            fine_smoother->set_R(R);  // Monotone

            auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>("mumps", "lu");
            // MultigridQR<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver, num_levels);  // QR
            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);  // Monotone

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
            interpolation_operators.resize(num_levels - 1);
            interpolation_operators[1] =
                std::make_shared<IPTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(QtIh));
            interpolation_operators[0] = std::make_shared<IPTransfer<Matrix, Vector>>(std::make_shared<Matrix>(Ih1));

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(verbose);
            // multigrid.verbose(true);

            // multigrid.mg_type(2);

            // This should be somewhere else...
            // multigrid.set_QR(Q, R);  // QR
            multigrid.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));

            multigrid.solve(QtAQ, Qtrhs, Qtx);
            x = Rot * Q * Qtx;
            // x = Rot*Qtx;

            // write("x.m", x);
            // write("IX.m", Qtx);
            // disp(x);
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(interior_point_qp_solver_test);
            UTOPIA_RUN_TEST(log_barrier_test);
            UTOPIA_RUN_TEST(log_barrier_qp_solver_test);
            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
            UTOPIA_RUN_TEST(MPRGP_test);
            UTOPIA_RUN_TEST(nblockgs_test);
        }

        void run_GS_QR() {
            if (mpi_world_size() > 1) return;

            print_backend_info();
            // UTOPIA_RUN_TEST(MG_QR_test);
        }

        QPSolverTest() : n(20) {}

        SizeType n = 20;
        bool verbose = false;
    };

    // FIXME merge with the other once it is poperly implemented
    template <class Matrix, class Vector>
    class PQPSolverTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;
        using IndexArray = typename Traits::IndexArray;
        using IndexSet = typename Traits::IndexSet;

        void run() {
            UTOPIA_RUN_TEST(MPRGP_DD);
            // FIXME
            UTOPIA_RUN_TEST(poly_qp);
        }

        void MPRGP_DD() {
            auto &&comm = Comm::get_default();

            // static const bool verbose = false;
            static const bool verbose = Traits::Backend == PETSC;

            Matrix A;
            Vector b;
            BoxConstraints<Vector> box;

            Vector oracle;
            std::stringstream c_ss;
            Chrono c;

            if (true) {
                c.start();

                SizeType n = 100;
                QPSolverTest<Matrix, Vector>::create_symm_lapl_test_data(comm, A, b, box, n, true);

                c.stop();
                c_ss << "Problem initialization\n" << c << "\n";

                // if (n <= 1e5) {
                //     c.start();
                //     Factorization<Matrix, Vector> solver;
                //     solver.solve(A, b, oracle);
                //     c.stop();
                //     c_ss << "Direct solver\n" << c << "\n";
                // }

            } else {
                c.start();

                // Path dir = "../data/test/CG_DD/mats_tests_2d_tri3";
                // Path dir = "../data/test/CG_DD/diffusion3d_P1_531k";
                Path dir = "../data/test/CG_DD/diffusion3d_P1_69k";
                // Path dir = "../data/test/CG_DD/diffusion2d_P2_103k";
                read(dir / "A", A);
                read(dir / "b", b);
                read(dir / "x", oracle);

                c.stop();
                c_ss << "Read problem from disk\n" << c << "\n";
            }

            ///////////////////////////////////////////////////////////////

            c.start();

            Vector x(layout(b), 0.);

            auto params = param_list(param("inner_solver",
                                           param_list(param("verbose", verbose),
                                                      param("atol", 1e-10),
                                                      param("rtol", 1e-10),
                                                      param("stol", 1e-10),
                                                      param("max_it", 2000))));

            // Test MFKSP
            if (false) {
                params.set("type", "ksp");
                params.set("pc_type", "none");
                params.set("ksp_type", "gmres");

                // params.set("type", "bcgs");
            }

            // params.set("use_preconditioner", true);
            // params.set("preconditioner_type", "amg");
            params.set("verbose", verbose);
            params.set("num_blocks", 3);

            // params.set("preconditioner_type", "inv");
            // params.set("use_preconditioner", false);

            ///////////////////////////////////////////////////////////////

            BDDLinearSolver<Matrix, Vector> solver;
            solver.read(params);
            solver.update(make_ref(A));

            c.stop();
            c_ss << "BDDLinearSolver::update\n" << c << "\n";

            ///////////////////////////////////////////////////////////////

            c.start();

            solver.apply(b, x);

            c.stop();
            c_ss << "BDDLinearSolver::solve\n" << c << "\n";

            ///////////////////////////////////////////////////////////////

            if (box.has_bound()) {
                BDDQPSolver<Matrix, Vector> qp_solver;

                /// Changing the input data
                box.lower_bound() = nullptr;
                auto &&u = *box.upper_bound();

                auto u_view = view_device(u);
                SizeType n_half = u.size() / 2;

                // Allow for unconstrained nodes
                parallel_for(
                    range_device(u), UTOPIA_LAMBDA(const SizeType i) {
                        if (i > n_half) {
                            u_view.set(i, 0.6);
                        }
                    });

                c.start();

                auto qp_params = param_list(param("infinity", 0.55),
                                            param("use_preconditioner", false),
                                            param("debug", verbose),
                                            param("inner_solver",
                                                  param_list(param("verbose", verbose),
                                                             param("atol", 1e-14),
                                                             param("rtol", 1e-14),
                                                             param("stol", 1e-14),
                                                             param("max_it", 2000))));

                qp_params.set("verbose", verbose);

                qp_solver.read(qp_params);

                qp_solver.set_box_constraints(box);
                qp_solver.update(make_ref(A));

                c.stop();
                c_ss << "BDDQPSolver::update\n" << c << "\n";

                ///////////////////////////////////////////////////////////////

                c.start();

                Vector x_qp(layout(b));
                qp_solver.apply(b, x_qp);

                c.stop();
                c_ss << "BDDLinearSolver::solve\n" << c << "\n";

                rename("x", x_qp);
                write("load_XQP.m", x_qp);
            }

            ///////////////////////////////////////////////////////////////

            if (verbose) {
                x.comm().root_print(c_ss.str());
            }

            ///////////////////////////////////////////////////////////////

            if (!empty(oracle)) {
                Scalar diff = norm2(x - oracle);

                if (diff > 1e-6 && verbose) {
                    comm.root_print(diff);

                    rename("o", oracle);
                    write("O.m", oracle);

                    rename("x", x);
                    write("X.m", x);
                }

                utopia_test_assert(diff < 1e-6);
            }
        }

        void poly_qp() {
            QPSolverTest<Matrix, Vector> s;

            OmniQPSolver<Matrix, Vector> solver;
            InputParameters in;

            in.set("backend", "any");
            in.set("type", "pg");
            // in.set("verbose", true);
            in.set("max_it", 2000);

            solver.read(in);

            s.run_qp_solver(solver);
        }
    };

    template <class Matrix, class Vector>
    class MonotoneMGTest {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        static void print_backend_info() {
            if (Utopia::instance().verbose() && mpi_world_rank() == 0) {
                utopia::out() << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        void monotone_amg_test() {
            // Does not work
            const static bool verbose = true;
            const static bool use_masks = false;
            int n_levels = 4;
            int n_coarse = 50;

            auto qp_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();
            auto agglomerator = std::make_shared<Agglomerate<Matrix>>();

            MonotoneAlgebraicMultigrid<Matrix, Vector> amg(qp_smoother, fine_smoother, direct_solver, agglomerator);
            amg.set_n_levels(n_levels);

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);
            auto funs = ml_problem.get_functions();

            Vector x, g;
            Matrix H;

            funs.back()->get_eq_constrains_values(x);
            funs.back()->gradient(x, g);
            funs.back()->hessian(x, H);

            Vector lower_bound(layout(g), -0.8), upper_bound(layout(g), 200.);

            amg.verbose(verbose);
            amg.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));

            utopia_test_assert(amg.solve(H, g, x));

            // rename("x", x);
            // write("X.m", x);
        }

        void monotone_mg_test() {
            const std::string data_path = Utopia::instance().get("data_path");

            const static bool verbose = false;
            const static bool use_masks = false;

            int n_levels = 6;
            int n_coarse = 50;

            using ProblemType = utopia::Poisson1D<Matrix, Vector>;
            MultiLevelTestProblem1D<Matrix, Vector, ProblemType> ml_problem(n_levels, n_coarse, !use_masks);
            auto funs = ml_problem.get_functions();

            Vector x, g;
            Matrix H;

            funs.back()->get_eq_constrains_values(x);
            funs.back()->gradient(x, g);
            funs.back()->hessian(x, H);

            auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            auto coarse_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            // auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            // auto coarse_smoother = std::make_shared<KSPSolver<Matrix, Vector>>();
            // coarse_smoother->pc_type("bjacobi");
            // coarse_smoother->ksp_type("cg");

            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();
            // auto direct_solver = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            MonotoneMultigrid<Matrix, Vector> multigrid(fine_smoother, coarse_smoother, direct_solver);

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
            interpolation_operators.resize(n_levels - 1);

            auto &transfers = ml_problem.get_transfer();
            for (SizeType i = 0; i < n_levels - 2; ++i) {
                interpolation_operators[i] = transfers[i];
            }

            auto t = std::static_pointer_cast<MatrixTransfer<Matrix, Vector>>(transfers[n_levels - 2]);
            interpolation_operators[n_levels - 2] = multigrid.new_fine_level_transfer(std::make_shared<Matrix>(t->I()));

            // std::make_shared<IPRTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(t->I()));

            Vector lower_bound(layout(g), -0.8), upper_bound(layout(g), 200.);

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(5);
            multigrid.post_smoothing_steps(5);
            multigrid.verbose(verbose);
            multigrid.set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));
            multigrid.update(make_ref(H));

            // avoids flip-floping of active nodes (and Galerkin assembly when nothing changes)
            multigrid.active_set().tol(1e-15);
            multigrid.apply(g, x);

            // disp(x);
        }

        void run() {
            print_backend_info();
            UTOPIA_RUN_TEST(monotone_mg_test);

            // FIXME
            // UTOPIA_RUN_TEST(monotone_amg_test);
        }
    };

    static void qp_solver() {
#ifdef UTOPIA_ENABLE_PETSC
        PQPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run_GS_QR();

        MonotoneMGTest<PetscMatrix, PetscVector>().run();
        // ProjectedGaussSeidelNewTest<PetscMatrix, PetscVector>().run();

#endif  // UTOPIA_ENABLE_PETSC

#ifdef UTOPIA_ENABLE_TRILINOS
        QPSolverTest<TpetraMatrixd, TpetraVectord>().run();
#endif  // UTOPIA_ENABLE_TRILINOS

#ifdef UTOPIA_ENABLE_BLAS
        QPSolverTest<BlasMatrixd, BlasVectord>().run();  // TODO(zulianp): : because blas is missing min operation ....
#endif                                                   // UTOPIA_ENABLE_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(qp_solver);
}  // namespace utopia
