#include "utopia_QPSolverTest.hpp"

#include "utopia.hpp"

#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGradient.hpp"

#include "test_problems/utopia_assemble_laplacian_1D.hpp"

#include "utopia_polymorphic_QPSolver.hpp"

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

        void MPRGP_test() const
        {
            std::cout<<"---- MPGRP solver --- \n"; 
            MPGRP<Matrix, Vector> qp_solver;
            run_qp_solver(qp_solver);


            SizeType n = 50; 

            Matrix A = sparse(n, n, 3);
            assemble_symmetric_laplacian_1D(A, true);

            auto h = 1./(n-1.);
            A = 1./h*A;

            {
                Range r = row_range(A);
                Write<Matrix> w(A);
                if(r.begin() == 0) {
                    A.set(0, 0, 1.);
                }

                if(r.end() == n) {
                    A.set(n-1, n-1, 1.);
                }
            }

            Vector  b = local_values(local_size(A).get(0), 50.0); 

            {
                Range row_range = range(b);
                Write<Vector> w(b);

                for(auto r = row_range.begin(); r!= row_range.end(); ++r)
                {
                    if(r >= n/2.)
                    {
                        b.set(r, -50.0); 
                    }
                    if(r == 0) 
                    {
                        b.set(r, 0); 
                    }

                    if(r == (n-1))
                    {
                        b.set(r, 0); 
                    }                    
                }
            }

            b = h*b; 

            Vector lb = local_values(local_size(A).get(0), -0.5); 
            Vector ub = local_values(local_size(A).get(0), 0.5); 

            Vector x = 0*b; 

            qp_solver.set_box_constraints(make_box_constaints(make_ref(lb), make_ref(ub)));
            qp_solver.verbose(verbose);
            qp_solver.max_it(100); 
            qp_solver.solve(A, b, x); 
        }


        void run()
        {
            print_backend_info();

            UTOPIA_RUN_TEST(pg_test);
            UTOPIA_RUN_TEST(pcg_test);
            UTOPIA_RUN_TEST(ngs_test);
            UTOPIA_RUN_TEST(MPRGP_test); 
        }

        QPSolverTest() : n(20) {}

        SizeType n = 20;
        bool verbose = false;
    };


    //FIXME merge with the other once it is poperly implemented
    template<class Matrix, class Vector>
    class PQPSolverTest {
    public:

        void run()
        {
           UTOPIA_RUN_TEST(poly_qp);
        }

        void poly_qp()
        {
            QPSolverTest<Matrix, Vector> s;

            PolymorphicQPSolver<Matrix, Vector> solver;
            InputParameters in;

            in.set("backend", "any");
            in.set("type",    "pg");
            // in.set("verbose", true);
            in.set("max-it",  2000);

            solver.read(in);

            s.run_qp_solver(solver);
        }

    };

    void run_qp_solver_test() {
        UTOPIA_UNIT_TEST_BEGIN("QPSolverTest");
#ifdef WITH_PETSC
        QPSolverTest<DSMatrixd, DVectord>().run();
        PQPSolverTest<DSMatrixd, DVectord>().run();

#endif //WITH_PETSC

#ifdef WITH_TRILINOS
        QPSolverTest<TSMatrixd, TVectord>().run();
#endif //WITH_TRILINOS

// #ifdef WITH_BLAS
//         QPSolverTest<Matrixd, Vectord>().run(); // TODO:: because blas is missing min operation .... 
// #endif //WITH_BLAS

        UTOPIA_UNIT_TEST_END("QPSolverTest");
    }
}
