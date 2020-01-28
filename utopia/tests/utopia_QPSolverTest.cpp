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
            MPGRP<Matrix, Vector> qp_solver;
            run_qp_solver(qp_solver);


            SizeType n = 100; 
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
            qp_solver.verbose(false);
            qp_solver.max_it(n*2); 
            qp_solver.set_eig_comp_tol(1e-1); 
            qp_solver.solve(A, b, x); 

            // disp(x, "x"); 
        }



        void ProjectedGS_QR()
        {

            std::cout<<"----------- Ciao Hardik ------------  \n"; 

            Vector rhs, x;
            Vector upper_bound,lower_bound;
            Matrix A, R, Q, Ih_fine;//, Ih1, Ih0; 

            const std::string data_path = Utopia::instance().get("data_path");

            read(data_path + "/forQR/b", rhs);
            read(data_path + "/forQR/x", x);
            read(data_path + "/forQR/A", A);            
            read(data_path + "/forQR/Q", Q);
            read(data_path + "/forQR/R", R);
            read(data_path + "/forQR/ub", upper_bound);
            read(data_path + "/forQR/lb", lower_bound);

            read(data_path + "/forQR/Ih", Ih_fine);
            //read(data_path + "/forQR/I2h", Ih1);
            //read(data_path + "/forQR/I3h", Ih0);


            auto num_levels = 2;


            // chop_abs(Q, 1e-4); 


            std::cout<<"A: "<< local_size(A).get(0) << "  \n"; 
            std::cout<<"Q: "<< local_size(Q).get(0) << "  \n"; 
            std::cout<<"R: "<< local_size(R).get(0) << ","<<local_size(R).get(1) <<   "  \n"; 
            std::cout<<"I: "<< local_size(Ih_fine).get(0) << "  \n"; 
            std::cout<<"rhs: "<< local_size(rhs).get(0) << "  \n"; 


            R = transpose(R);   


            Matrix QtAQ  = transpose(Q)*A*Q;
            Matrix QtIh  = transpose(Q)*Ih_fine;
            Vector Qtrhs = transpose(Q)*rhs;
            Vector Qtx   = transpose(Q)*x; 


            // disp(rhs);
            // disp(A);
            
            // Vector x = local_values(local_size(rhs).get(0), 0.0);


            auto smoother_fine = std::make_shared<ProjectedGaussSeidelQR<Matrix, Vector>>();
            auto direct_solver = std::make_shared<Factorization<Matrix, Vector> >("mumps", "lu");

            //Vector upper_bound = local_values(local_size(rhs).get(0), 0);
            //Vector lower_bound  = local_values(local_size(rhs).get(0), 0);

            // smoother->max_it(20);
            // smoother->stol(1e-14);
            // smoother->verbose(true);
            smoother_fine->max_it(100);
            smoother_fine->n_local_sweeps(1);            
            smoother_fine->set_box_constraints(make_box_constaints(make_ref(lower_bound),  make_ref(upper_bound)));
            smoother_fine->set_R(R);
            smoother_fine->verbose(true);
            // smoother_fine->solve(QtAQ, Qtrhs, Qtx);
            // x = Q * Qtx; 

            //write("x.m", x); 
            //write("IX.m", Qtx);
            // exit(0);

            // MG test starts here...
            // std::vector<std::shared_ptr <Matrix> > interpolation_operators;
            // interpolation_operators.push_back(make_ref(QtIh));

            std::vector<std::shared_ptr<Transfer<Matrix, Vector> > > interpolation_operators;            
            // interpolation_operators.resize(3);
            // interpolation_operators[2] = std::make_shared<MatrixTruncatedTransfer<Matrix, Vector> >(std::make_shared<Matrix>(QtIh));
            // interpolation_operators[1] = std::make_shared<MatrixTruncatedTransfer<Matrix, Vector> >(std::make_shared<Matrix>(Ih1));
            // interpolation_operators[0] = std::make_shared<MatrixTruncatedTransfer<Matrix, Vector> >(std::make_shared<Matrix>(Ih0));
    
            interpolation_operators.resize(1);
            interpolation_operators[0] = std::make_shared<MatrixTruncatedTransfer<Matrix, Vector> >(std::make_shared<Matrix>(QtIh));
            

            auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();
            MultigridQR<Matrix, Vector> multigrid(smoother_fine, direct_solver, num_levels);

            //multigrid.set_smoother(smoother_fine, num_levels-1); 
            multigrid.set_transfer_operators(interpolation_operators);
            multigrid.fix_semidefinite_operators(true); 
            multigrid.max_it(10);
            multigrid.use_line_search(true); 
            multigrid.update(make_ref(QtAQ));
            multigrid.verbose(true);
            multigrid.pre_smoothing_steps(1);


            // This should be somewhere else... 
            multigrid.set_R(R); 
            multigrid.set_upper_bound(upper_bound); 
            multigrid.set_lower_bound(lower_bound); 



            multigrid.apply(Qtrhs, Qtx);
            x = Q * Qtx; 

            write("x.m", x); 
            write("IX.m", Qtx);
            // disp(x);

        }






        void run()
        {
            print_backend_info();

            // UTOPIA_RUN_TEST(pg_test);
            // UTOPIA_RUN_TEST(pcg_test);
            // UTOPIA_RUN_TEST(ngs_test);
            // UTOPIA_RUN_TEST(MPRGP_test); 


            UTOPIA_RUN_TEST(ProjectedGS_QR); 

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
