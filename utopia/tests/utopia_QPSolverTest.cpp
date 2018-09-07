#include "utopia_QPSolverTest.hpp"

#include "utopia.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "test_problems/utopia_assemble_laplacian_1D.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class QPSolverTest {
    public:

        static void print_backend_info()
        {
            if(Utopia::instance().verbose() && mpi_world_rank() == 0) {
                std::cout << "\nBackend: " << backend_info(Vector()).get_name() << std::endl;
            }
        }

        template<class QPSolver>
        void run_qp_solver(QPSolver &qp_solver) const {

            Matrix m = sparse(n, n, 3);
            assemble_laplacian_1D(m);
            {
                Range r = row_range(m);
                Write<Matrix> w(m);
                if(r.begin() == 0) {
                    m.set(0, 0, 1.);
                    m.set(0, 1, 0);
                }

                if(r.end() == n) {
                    m.set(n-1, n-1, 1.);
                    m.set(n-1, n-2, 0);
                }
            }

            Vector rhs = values(n, 1.);
            {
                //Creating test vector (alternative way see [assemble vector alternative], which might be easier for beginners)
                Range r = range(rhs);
                Write<Vector> w(rhs);

                if(r.begin() == 0) {
                    rhs.set(0, 0);
                }

                if(r.end() == n) {
                    rhs.set(n-1, 0.);
                }
            }

            Vector upper_bound = values(n, 100.0);
            Vector solution    = zeros(n);


            qp_solver.max_it(n*40);
            qp_solver.verbose(verbose);
            qp_solver.set_box_constraints(make_upper_bound_constraints(make_ref(upper_bound)));

            Chrono c;
            c.start();
            bool ok = qp_solver.solve(m, rhs, solution);
            c.stop();

            utopia_test_assert(ok);
        }

        void pg_test() const
        {
            ProjectedGradient<Matrix, Vector> pg;
            run_qp_solver(pg);
        }

        void pcg_test() const
        {
            ProjectedConjugateGradient<Matrix, Vector> pcg;
            run_qp_solver(pcg);
        }

        void ngs_test() const
        {
            ProjectedGaussSeidel<Matrix, Vector> pgs;
            run_qp_solver(pgs);
        }

        void run()
        {
            print_backend_info();

            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
        }

        SizeType n = 20;
        bool verbose = false;
    };

    void run_qp_solver_test() {
        UTOPIA_UNIT_TEST_BEGIN("QPSolverTest");
#ifdef WITH_PETSC
        QPSolverTest<DSMatrixd, DVectord>().run();
#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        QPSolverTest<TSMatrixd, TVectord>().run();
#endif //WITH_TRILINOS

#ifdef WITH_BLAS
        QPSolverTest<Matrixd, Vectord>().run();
#endif //WITH_BLAS

        UTOPIA_UNIT_TEST_END("QPSolverTest");
    }
}
