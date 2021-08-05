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

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"
#endif  // UTOPIA_WITH_PETSC

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolverTest {
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

        static void create_symm_lapl_test_data(Comm &comm, Matrix &A, Vector &b, BoxConstraints<Vector> &box) {
            SizeType n = 100;
            A.sparse(layout(comm, Traits::decide(), Traits::decide(), n, n), 3, 2);
            assemble_symmetric_laplacian_1D(A, true);

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

            b *= 0.5;
            // box.lower_bound() = nullptr;
            // box.upper_bound() = nullptr;
            QuadraticFunction<Matrix, Vector> fun(make_ref(A), make_ref(b));

            InputParameters params;
            // params.set("verbose", true);
            params.set("barrier_parameter", 1e-5);
            params.set("barrier_parameter_shrinking_factor", 0.7);

            LogBarrierFunction<Matrix, Vector> barrier(make_ref(fun), make_ref(box));
            barrier.read(params);

            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            cg.set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector>>());
            // cg.verbose(true);
            cg.max_it(20);

            Newton<Matrix, Vector> newton(make_ref(cg));
            // newton.verbose(true);

            Vector x(layout(b));
            newton.solve(barrier, x);

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

            UTOPIA_RUN_TEST(log_barrier_test);
            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
            UTOPIA_RUN_TEST(MPRGP_test);
            UTOPIA_RUN_TEST(nblockgs_test);
        }

        void run_GS_QR() {
            if (mpi_world_size() > 1) return;

            print_backend_info();
            UTOPIA_RUN_TEST(MG_QR_test);
        }

        QPSolverTest() : n(20) {}

        SizeType n = 20;
        bool verbose = false;
    };

    // FIXME merge with the other once it is poperly implemented
    template <class Matrix, class Vector>
    class PQPSolverTest {
    public:
        void run() {
            // FIXME
            UTOPIA_RUN_TEST(poly_qp);
        }

        void poly_qp() {
            QPSolverTest<Matrix, Vector> s;

            OmniQPSolver<Matrix, Vector> solver;
            InputParameters in;

            in.set("backend", "any");
            in.set("type", "pg");
            // in.set("verbose", true);
            in.set("max-it", 2000);

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

            rename("x", x);
            write("X.m", x);
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
#ifdef UTOPIA_WITH_PETSC
        QPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run_GS_QR();

        PQPSolverTest<PetscMatrix, PetscVector>().run();
        MonotoneMGTest<PetscMatrix, PetscVector>().run();
        // ProjectedGaussSeidelNewTest<PetscMatrix, PetscVector>().run();

#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
        QPSolverTest<TpetraMatrixd, TpetraVectord>().run();
#endif  // UTOPIA_WITH_TRILINOS

#ifdef UTOPIA_WITH_BLAS
        QPSolverTest<BlasMatrixd, BlasVectord>().run();  // TODO(zulianp): : because blas is missing min operation ....
#endif                                                   // UTOPIA_WITH_BLAS
    }

    UTOPIA_REGISTER_TEST_FUNCTION(qp_solver);
}  // namespace utopia
