#include "utopia_Testing.hpp"

#include "utopia.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "test_problems/utopia_QPSolverTestProblem.hpp"
#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_ProjectedGaussSeidelNew.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

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

        void MPRGP_test() const {
            MPGRP<Matrix, Vector> qp_solver;
            run_qp_solver(qp_solver);

            auto &&comm = Comm::get_default();

            SizeType n = 100;
            Matrix A;
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

            Vector b(row_layout(A), 50.0);

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

            Vector lb(row_layout(A), -0.5);
            Vector ub(row_layout(A), 0.5);

            Vector x = 0 * b;

            qp_solver.set_box_constraints(make_box_constaints(make_ref(lb), make_ref(ub)));
            qp_solver.verbose(false);
            qp_solver.max_it(n * 2);
            qp_solver.set_eig_comp_tol(1e-1);
            qp_solver.solve(A, b, x);
        }

        void ProjectedGS_QR() {
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

            auto num_levels = 3;

            // chop_abs(Q, 1e-7);

            // std::cout<<"A: "<< local_size(A).get(0) << "  \n";
            // std::cout<<"Q: "<< local_size(Q).get(0) << "  \n";
            // std::cout<<"R: "<< local_size(R).get(0) << ","<<local_size(R).get(1) <<   "  \n";
            // std::cout<<"I: "<< local_size(Ih_fine).get(0) << "  \n";
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

            auto smoother_fine = std::make_shared<ProjectedGaussSeidelQR<Matrix, Vector>>();
            auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>("mumps", "lu");
            MultigridQR<Matrix, Vector> multigrid(smoother_fine, coarse_smoother, direct_solver, num_levels);

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> interpolation_operators;
            interpolation_operators.resize(num_levels - 1);
            interpolation_operators[1] =
                std::make_shared<MatrixTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(QtIh));
            interpolation_operators[0] =
                std::make_shared<MatrixTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(Ih1));

            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.max_it(40);
            multigrid.pre_smoothing_steps(3);
            multigrid.post_smoothing_steps(3);
            multigrid.verbose(false);
            // multigrid.mg_type(2);

            // This should be somewhere else...
            multigrid.set_QR(Q, R);
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

            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
            UTOPIA_RUN_TEST(MPRGP_test);
        }

        void run_GS_QR() {
            if (mpi_world_size() > 1) return;

            print_backend_info();
            UTOPIA_RUN_TEST(ProjectedGS_QR);
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

            PolymorphicQPSolver<Matrix, Vector> solver;
            InputParameters in;

            in.set("backend", "any");
            in.set("type", "pg");
            // in.set("verbose", true);
            in.set("max-it", 2000);

            solver.read(in);

            s.run_qp_solver(solver);
        }
    };

    // template<class Matrix, class Vector>
    // class ProjectedGaussSeidelNewTest {
    // public:
    //     void run()
    //     {
    //         //FIXME
    //        UTOPIA_RUN_TEST(pg_new_test);
    //     }

    //     template<class QPSolver>
    //     void run_qp_solver(QPSolver &qp_solver) const {
    //         QPSolverTestProblem<Matrix, Vector>::run(n, verbose, qp_solver, true);
    //     }

    //     void pg_new_test()
    //     {
    //         ProjectedGaussSeidelNew<Matrix, Vector> pgs;
    //         run_qp_solver(pgs);
    //     }

    //     SizeType n = 20;
    //     bool verbose = true;
    // };

    static void qp_solver() {
#ifdef UTOPIA_WITH_PETSC
        QPSolverTest<PetscMatrix, PetscVector>().run();
        QPSolverTest<PetscMatrix, PetscVector>().run_GS_QR();

        PQPSolverTest<PetscMatrix, PetscVector>().run();
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
