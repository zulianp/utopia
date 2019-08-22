#include "utopia.hpp"
#include "test_problems/utopia_TestProblems.hpp"

namespace utopia
{

    template<typename Matrix, typename Vector>
    class HckTests 
    {

    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        HckTests(const SizeType & n, const SizeType & n_levels=2, const Scalar & lambda = 1.0, const bool verbose = false, const bool output_flg = false): 
        n_(n), 
        n_levels_(n_levels),
        lambda_(lambda), 
        verbose_(verbose),
        output_vtk_(output_flg)
        {
            input_params_.set("atol", 1e-7);
            input_params_.set("rtol", 1e-10);
            input_params_.set("stol", 1e-10);        
            input_params_.set("verbose", verbose_); 
        }

        void run()
        {
            // UTOPIA_RUN_TEST(newton_test);            

            // UTOPIA_RUN_TEST(STCG_test); 
            // UTOPIA_RUN_TEST(MPGRP); 

            // UTOPIA_RUN_TEST(TR_unconstrained);
            // UTOPIA_RUN_TEST(TR_constrained); 

            // UTOPIA_RUN_TEST(QuasiTR_unconstrained);
            // UTOPIA_RUN_TEST(QuasiTR_constrained);


            UTOPIA_RUN_TEST(RMTR_unconstrained); 
        }

        template<class QPSolverTemp>
        void QP_solve(QPSolverTemp &qp_solver) const 
        {
            Bratu2D<Matrix, Vector> fun(n_);
            Vector x = fun.initial_guess(); 
            Matrix H;
            Vector g; 

            fun.hessian(x, H); 
            fun.gradient(x, g); 
            x *= 0.0;


            if(dynamic_cast<QPSolver<Matrix, Vector> *>(qp_solver.get()))
            {
                QPSolver<Matrix, Vector> * qp_box = dynamic_cast<QPSolver<Matrix, Vector> *>(qp_solver.get());
                Vector lb = local_values(local_size(x).get(0), -9e9); 
                Vector ub = local_values(local_size(x).get(0), 9e9); 
                qp_box->set_box_constraints(make_box_constaints(make_ref(lb), make_ref(ub)));      
                qp_box->solve(H, -1.0*g, x); 
            }
            else
            {
                qp_solver->solve(H, -1.0*g, x); 
            }
        }


        void STCG_test()
        {
            auto QP_solver = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            QP_solver->set_preconditioner(std::make_shared<InvDiagPreconditioner<Matrix, Vector> >());
            QP_solver->atol(1e-10);
            QP_solver->max_it(n_*n_);
            QP_solver->verbose(true); 
            QP_solver->current_radius(9e9); 

            QP_solve(QP_solver); 
        }

        void MPGRP()
        {
            auto QP_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            QP_solver->atol(1e-10);
            QP_solver->max_it(n_*n_);
            QP_solver->verbose(false); 
    
            QP_solve(QP_solver); 
        }


        void TR_unconstrained()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            fun.describe(); 

            auto subproblem = std::make_shared<utopia::Lanczos<Matrix, Vector> >();
            // auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            subproblem->pc_type("bjacobi"); 
            subproblem->atol(1e-14);
            subproblem->max_it(100000);
            
            TrustRegion<Matrix, Vector> tr_solver(subproblem);
            tr_solver.read(input_params_); 
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }

        void newton_test()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            fun.describe(); 

            auto lsolver = std::make_shared<GMRES<Matrix, Vector> >();
            lsolver->pc_type("bjacobi"); 
            
            Newton<Matrix, Vector> solver(lsolver);
            solver.read(input_params_); 
            solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }

        void QuasiTR_unconstrained()
        {
            Bratu2D<Matrix, Vector> fun(n_);
            Vector x = fun.initial_guess(); 
            SizeType memory_size = 5; 

            auto subproblem = std::make_shared<utopia::SteihaugToint<Matrix, Vector, HOMEMADE> >();
            subproblem->set_preconditioner(std::make_shared<IdentityPreconditioner<Vector> >());
            subproblem->atol(1e-14);
            subproblem->max_it(100000);

            auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);   
            hess_approx->theta_min(1.0);
            hess_approx->damping_tech(POWEL); 
            hess_approx->scaling_tech(ADAPTIVE); 

            QuasiTrustRegion<Vector> tr_solver(hess_approx, subproblem);
            tr_solver.read(input_params_); 
            tr_solver.max_it(5000);
            tr_solver.delta0(10000);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x, "QuasiTRUnstrained.vtk");
        }        


        void TR_constrained()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            fun.describe(); 

            auto qp_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            qp_solver->atol(1e-10);
            qp_solver->max_it(n_*n_);
            qp_solver->verbose(false); 

            // auto lsolver = std::make_shared<GMRES<Matrix, Vector> >();
            // auto qp_solver =  std::make_shared<utopia::TaoQPSolver<Matrix, Vector> >(lsolver);
            // qp_solver->atol(1e-11);            
            
            TrustRegionVariableBound<Matrix, Vector> tr_solver(qp_solver);
            
            Vector ub, lb; 
            fun.upper_bound(ub); 
            fun.lower_bound(lb);             
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));

            tr_solver.set_box_constraints(box);
            tr_solver.read(input_params_); 
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x);
        }        



        void QuasiTR_constrained()
        {
            Bratu2D<Matrix, Vector> fun(n_);
            Vector x = fun.initial_guess(); 
            SizeType memory_size = 10; 

            auto qp_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            qp_solver->atol(1e-10);
            qp_solver->max_it(n_*n_);

            auto hess_approx   = std::make_shared<LBFGS<Vector> >(memory_size);   
            // hess_approx->theta_min(1.0);
            // hess_approx->damping_tech(POWEL); 
            // hess_approx->scaling_tech(ADAPTIVE); 

            Vector ub, lb; 
            fun.upper_bound(ub); 
            fun.lower_bound(lb);             
            auto box = make_box_constaints(make_ref(lb), make_ref(ub));       


            QuasiTrustRegionVariableBound<Vector>  tr_solver(hess_approx, qp_solver);
            tr_solver.set_box_constraints(box);
            tr_solver.read(input_params_); 
            tr_solver.delta0(0.01);
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x, "QuasiTRConstrained.vtk");
        }    


        void RMTR_unconstrained()
        {
            // auto n_coarse = 50; 
            auto n_coarse = n_; 

            Petsc2DMultilevelTestProblem<Matrix, Vector, Bratu2D<Matrix, Vector> > multilevel_problem(n_levels_, n_coarse); 

            auto fun = multilevel_problem.level_functions_[n_levels_-1];
            Bratu2D<Matrix, Vector> * fun_Bratu2D = dynamic_cast<Bratu2D<Matrix, Vector> *>(fun.get());
            Vector x = fun_Bratu2D->initial_guess(); 
            fun_Bratu2D->describe(); 

            auto tr_strategy_coarse = std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            tr_strategy_coarse->pc_type("lu"); 
            tr_strategy_coarse->atol(1e-14);

            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >();
            tr_strategy_fine->pc_type("jacobi"); 
            tr_strategy_fine->atol(1e-14);

            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, FIRST_ORDER> >(n_levels_);
            auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(n_levels_);

            
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);                        

            rmtr->set_transfer_operators(multilevel_problem.transfers_);

            rmtr->max_it(10);
            rmtr->max_coarse_it(3);
            rmtr->max_QP_coarse_it(300);

            rmtr->max_smoothing_it(2);
            rmtr->max_QP_smoothing_it(20); 


            rmtr->delta0(1e9);
            rmtr->atol(1e-6);
            rmtr->rtol(1e-10);
            rmtr->set_grad_smoothess_termination(0.000001);
            rmtr->set_eps_grad_termination(1e-7);

            rmtr->verbose(verbose_);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            rmtr->set_functions( multilevel_problem.level_functions_);
            rmtr->handle_equality_constraints();            
            rmtr->solve(x);
        }


    private: 
        SizeType n_; 
        SizeType n_levels_; 
        Scalar lambda_; 
        bool verbose_; 
        bool output_vtk_; 

        InputParameters input_params_;

    };

    void runHckTest()
    {
        UTOPIA_UNIT_TEST_BEGIN("HckTests");
        #ifdef  WITH_PETSC

            auto n_levels = 4; 
            auto coarse_dofs = 125; 

            HckTests<DSMatrixd, DVectord>(coarse_dofs, n_levels, 1.0, true, true).run();
        #endif
        UTOPIA_UNIT_TEST_END("HckTests");

    }
}
