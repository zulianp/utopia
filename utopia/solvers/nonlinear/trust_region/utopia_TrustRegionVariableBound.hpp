#ifndef UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP
#define UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP

#include "utopia_NewtonBase.hpp"
#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_VariableBoundSolverInterface.hpp"

 namespace utopia 
 {
    	template<class Matrix, class Vector>
      /**
       * @brief      Trust region solver taking into account also bound constraints.
       */ 
     	class TrustRegionVariableBound :  public VariableBoundSolverInterface<Vector>, 
                                        public TrustRegionBase<Matrix, Vector>, 
                                        public NewtonBase<Matrix, Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRBoxSubproblem<Matrix, Vector>       TRBoxSubproblem; 
        typedef utopia::TrustRegionBase<Matrix, Vector>       TrustRegionBase; 
        typedef utopia::NewtonBase<Matrix, Vector>       NonLinearSolver;
     	
     	public:                                                                       // once generic, then = std::shared_ptr<ProjectedGaussSeidel<Matrix, Vector> >()
      TrustRegionVariableBound( const std::shared_ptr<TRBoxSubproblem> &tr_subproblem,
                                const Parameters params = Parameters()) : 
                                NonLinearSolver(tr_subproblem, params), it_successful_(0)  
      {
        set_parameters(params);        
      }

      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
      virtual void set_parameters(const Parameters params) override
      {

        NonLinearSolver::set_parameters(params);
        TrustRegionBase::set_parameters(params);
      }

      /**
       * @brief      Trust region solve. 
       *
       * @param      fun   The nonlinear solve function.
       * @param      x_k   Initial gues/ solution
       *
       *  
       * @return     true
       */
      bool solve(Function<Matrix, Vector> &fun, Vector &x_k) override 
      {
        using namespace utopia;

        // passing solver and parameters into subproblem 
        bool converged = false, accepted = true; 
        NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

        Scalar delta, ared, pred, rho, E_old, E_new; 

        SizeType it = 0;
        it_successful_ = 0; 
        const Scalar infty = std::numeric_limits<Scalar>::infinity();
        Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;
        bool rad_flg = false; 

        Vector g = 0*x_k, p_k = 0*x_k, x_k1 = 0*x_k;
        Matrix H; 

        fun.gradient(x_k, g);

        this->make_iterate_feasible(x_k); 

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg); 
        // delta = 10;        // testing 

        g0_norm = norm2(g);
        g_norm = g0_norm;
        
        // print out - just to have idea how we are starting 
        if(this->verbose_)
        {
          this->init_solver("TRUST_REGION_BASE",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "ared","pred", "rho", "delta_k", "|| p_k || "});
          PrintInfo::print_iter_status(it, {g_norm}); 
        }

        it++; 
        fun.value(x_k, E_old); 
        fun.gradient(x_k, g);

        // solve starts here 
        while(!converged)
        {
          if(accepted)
          {
            fun.value(x_k, E_old); 
            fun.hessian(x_k, H); 
          }

    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------          
          if(TRBoxSubproblem * tr_subproblem = dynamic_cast<TRBoxSubproblem*>(this->linear_solver_.get()))
          {
            p_k = 0 * p_k; 
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta);
            tr_subproblem->tr_constrained_solve(H, g, p_k, box);
          }


          pred = this->get_pred(g, H, p_k); 
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
          // trial point 
          x_k1 = x_k + p_k; 
          // this->make_iterate_feasible(x_k1); 


          // value of the objective function with correction 
          fun.value(x_k1, E_new);

          // decrease ratio 
          ared = E_old - E_new;                // reduction observed on objective function
          pred = std::abs(pred); 
          rho = ared/ pred;               // decrease ratio         

    //----------------------------------------------------------------------------
    //     acceptance of trial point 
    //----------------------------------------------------------------------------
          if(ared < 0 || pred < 0)
            rho = 0.0; 
          else if(rho != rho)
            rho = 0.0; 


          this->trial_point_acceptance(rho, x_k1, x_k); 
          
          if (rho >= this->rho_tol())
            it_successful_++; 

    //----------------------------------------------------------------------------
    //    convergence check 
    //----------------------------------------------------------------------------
          if(accepted)
            fun.gradient(x_k, g); 
          
          g_norm = this->criticality_measure_infty(x_k, g); 
          s_norm = norm2(p_k); 

          if(this->verbose_)
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, ared, pred, rho, delta, s_norm}); 

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta); 

    //----------------------------------------------------------------------------
    //      tr. radius update 
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta, true); 
          it++; 
        }

        // some benchmarking 
        this->print_statistics(it);      

          return false;
      }


      virtual void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &ls) override
      {
          auto linear_solver = this->linear_solver(); 
          if (dynamic_cast<TRBoxSubproblem *>(linear_solver.get()) != nullptr)
          {
              TRBoxSubproblem * tr_sub = dynamic_cast<TRBoxSubproblem *>(linear_solver.get());
              tr_sub->set_linear_solver(ls);
          }
      }

      virtual void set_trust_region_strategy(const std::shared_ptr<TRBoxSubproblem> &tr_linear_solver)
      {
        NonLinearSolver::set_linear_solver(tr_linear_solver); 
      }


    protected: 

        virtual void print_statistics(const SizeType & it) override
        {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if(!non_data_path.empty())
            {
                CSVWriter writer;
                if (mpi_world_rank() == 0)
                {
                    if(!writer.file_exists(non_data_path))
                    {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"num_its", "it_successful", "time"});
                    }
                    else
                        writer.open_file(non_data_path);
                    
                    writer.write_table_row<Scalar>({Scalar(it), Scalar(it_successful_),  this->get_time()});
                    writer.close_file();
                }
            }
        }


    private:
      SizeType it_successful_; 


  };

}

#endif //UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP

