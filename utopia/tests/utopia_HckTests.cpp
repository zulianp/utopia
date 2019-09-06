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
            input_params_.set("max-it", 10); 

            // RMTR specific parameters 
            input_params_.set("max_coarse_it", 2); 
            input_params_.set("max_sucessful_coarse_it", 1); 
            input_params_.set("max_QP_coarse_it", 300);             
            input_params_.set("pre_smoothing_steps", 10); 
            input_params_.set("post_smoothing_steps", 10); 
            input_params_.set("max_sucessful_smoothing_it", 1);   
            input_params_.set("max_QP_smoothing_it", 1);   
            input_params_.set("delta0", 1);   
            input_params_.set("grad_smoothess_termination", 1e-8);      
        }

        void run()
        {
            UTOPIA_RUN_TEST(newton_test);            

            // UTOPIA_RUN_TEST(STCG_test); 
            // UTOPIA_RUN_TEST(MPGRP); 

            UTOPIA_RUN_TEST(TR_unconstrained);
            UTOPIA_RUN_TEST(TR_constrained); 
            
            UTOPIA_RUN_TEST(Poisson_test); 

            UTOPIA_RUN_TEST(QuasiTR_unconstrained);
            UTOPIA_RUN_TEST(QuasiTR_constrained);


            UTOPIA_RUN_TEST(RMTR_unconstrained); 
            UTOPIA_RUN_TEST(RMTR_l2_linear); 
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

        void Poisson_test()
        {
            Poisson3D<Matrix, Vector> fun(10);

            Vector b; 
            Matrix H; 
            Vector x = fun.initial_guess(); 
            
            if(verbose_)
                fun.describe(); 

            auto subproblem = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "sor", false);
            subproblem->atol(1e-14);
            subproblem->max_it(1000);
            
            TrustRegion<Matrix, Vector> tr_solver(subproblem);
            tr_solver.read(input_params_); 
            tr_solver.delta0(0.01); 
            tr_solver.atol(1e-12); 
            tr_solver.solve(fun, x);


            PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > multilevel_problem(3, n_levels_, n_); 

            auto fun1 = multilevel_problem.level_functions_[n_levels_-1];
            Poisson3D<Matrix, Vector> * fun_Laplace = dynamic_cast<Poisson3D<Matrix, Vector> *>(fun1.get());


            x = fun_Laplace->initial_guess(); 
            fun_Laplace->gradient(x, b); 
            fun_Laplace->hessian(x, H); 


            auto direct_solver = std::make_shared<utopia::RedundantLinearSolver<Matrix, Vector> >();
            direct_solver->number_of_parallel_solves(mpi_world_size()); 

            auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
            
            Multigrid<Matrix, Vector> multigrid(smoother, direct_solver);
            multigrid.set_transfer_operators(multilevel_problem.transfers_);
            multigrid.update(make_ref(H));
            multigrid.read(input_params_);
            multigrid.apply(b, x);            


            if(output_vtk_)
                fun.output_to_VTK(x);
        }


        void TR_unconstrained()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            
            if(verbose_)
                fun.describe(); 

            auto subproblem = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
            
            TrustRegion<Matrix, Vector> tr_solver(subproblem);
            tr_solver.read(input_params_); 

            if(output_vtk_)
                fun.output_to_VTK(x);
        }

        void newton_test()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            
            if(verbose_)
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
            tr_solver.solve(fun, x);

            if(output_vtk_)
                fun.output_to_VTK(x, "QuasiTRUnstrained.vtk");
        }        


        void TR_constrained()
        {
            Bratu2D<Matrix, Vector> fun(n_, 5.0);
            Vector x = fun.initial_guess(); 
            
            if(verbose_)
                fun.describe(); 

            auto qp_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
            qp_solver->atol(1e-10);
            qp_solver->max_it(n_*n_);
            qp_solver->verbose(false); 
       

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
            PetscMultilevelTestProblem<Matrix, Vector, Bratu2D<Matrix, Vector> > multilevel_problem(2, n_levels_, n_); 

            auto fun = multilevel_problem.level_functions_[n_levels_-1];
            Bratu2D<Matrix, Vector> * fun_Bratu2D = dynamic_cast<Bratu2D<Matrix, Vector> *>(fun.get());
            Vector x = fun_Bratu2D->initial_guess(); 

            if(verbose_)
                fun_Bratu2D->describe(); 

            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");


            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER> >(n_levels_);
            auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(n_levels_);

            // Set TR-QP strategies 
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);                        

            // Transfers and objective functions
            rmtr->set_transfer_operators(multilevel_problem.transfers_);
            rmtr->set_functions( multilevel_problem.level_functions_);    

 
            rmtr->norm_schedule(NormSchedule::OUTER_CYCLE);
            rmtr->read(input_params_); 
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
                
            // Solve 
            rmtr->solve(x);


            if(output_vtk_)
                fun_Bratu2D->output_to_VTK(x, "RMTR_output.vtk");            
        }

        void RMTR_l2_linear()
        {

            PetscMultilevelTestProblem<Matrix, Vector, Poisson3D<Matrix, Vector> > multilevel_problem(3, n_levels_, n_); 

            auto fun = multilevel_problem.level_functions_[n_levels_-1];
            Poisson3D<Matrix, Vector> * fun_Poisson3D = dynamic_cast<Poisson3D<Matrix, Vector> *>(fun.get());
            Vector x = fun_Poisson3D->initial_guess(); 

            if(verbose_)
                fun_Poisson3D->describe(); 


            auto tr_strategy_coarse = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg", "lu", true);
            auto tr_strategy_fine = std::make_shared<utopia::Lanczos<Matrix, Vector> >("sor");
            

            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, GALERKIN> >(n_levels_);
            // auto rmtr = std::make_shared<RMTR<Matrix, Vector, SECOND_ORDER> >(n_levels_);
            auto rmtr = std::make_shared<RMTR<Matrix, Vector, FIRST_ORDER> >(n_levels_);

            // Set TR-QP strategies 
            rmtr->set_coarse_tr_strategy(tr_strategy_coarse);
            rmtr->set_fine_tr_strategy(tr_strategy_fine);                        

            // Transfers and objective functions
            rmtr->set_transfer_operators(multilevel_problem.transfers_);
            rmtr->set_functions( multilevel_problem.level_functions_);    


            rmtr->norm_schedule(NormSchedule::OUTER_CYCLE);
            // rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
            rmtr->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);
            rmtr->read(input_params_); 
                
            // Solve 
            rmtr->solve(x);


            if(output_vtk_)
                fun_Poisson3D->output_to_VTK(x, "RMTR__linear_output.vtk");            
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

            auto n_levels = 3; 
            auto coarse_dofs = 5; 

            // auto n_levels = 4; 
            // auto coarse_dofs = 100; 

            HckTests<DSMatrixd, DVectord>(coarse_dofs, n_levels, 1.0, false, false).run();
        #endif
        UTOPIA_UNIT_TEST_END("HckTests");

    }
}
