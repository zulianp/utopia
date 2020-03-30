#ifndef UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP
#define UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

#include "utopia_NonLinearSolver.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"

 namespace utopia
 {
    	template<class Vector>
     	class QuasiTrustRegionVariableBound final:  public VariableBoundSolverInterface<Vector>,
                                                  public TrustRegionBase<Vector>,
                                                  public QuasiNewtonBase<Vector>
      {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        typedef utopia::MatrixFreeQPSolver<Vector>    MatrixFreeQPSolver;
        typedef utopia::TrustRegionBase<Vector>       TrustRegionBase;
        typedef utopia::QuasiNewtonBase<Vector>       NonLinearSolver;

        typedef utopia::HessianApproximation<Vector>  HessianApproximation;

    public:
      QuasiTrustRegionVariableBound(const std::shared_ptr <HessianApproximation> &hessian_approx,
                                    const std::shared_ptr <MatrixFreeQPSolver> &tr_subproblem) :
                                    NonLinearSolver(hessian_approx, tr_subproblem),
                                    it_successful_(0),
                                    initialized_(false),
                                    layout_(0)

      {

      }

      void read(Input &in) override
      {
        TrustRegionBase::read(in);
        NonLinearSolver::read(in);
      }

      void print_usage(std::ostream &os) const override
      {
        TrustRegionBase::print_usage(os);
        NonLinearSolver::print_usage(os);
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
      bool solve(FunctionBase<Vector> &fun, Vector &x_k) override
      {
        using namespace utopia;

        // passing solver and parameters into subproblem
        bool converged = false;
        NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

        Scalar delta, ared, pred, rho, E_old, E_new;

        this->fill_empty_bounds(local_size(x_k));
        this->make_iterate_feasible(x_k);

        SizeType it = 0;
        it_successful_ = 0;
        const Scalar infty = std::numeric_limits<Scalar>::infinity();
        Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;
        bool rad_flg = false;


        SizeType layout_x = local_size(x_k);
        if(!initialized_ || !g.comm().conjunction(layout_ == layout_x))
        {
          init_memory(layout_x);
        }

        fun.gradient(x_k, g);

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);

        // this->initialize_approximation(x_k, g);
        QuasiNewtonBase<Vector>::init_memory(x_k, g);
        auto multiplication_action = this->hessian_approx_strategy_->build_apply_H();


        g0_norm = this->criticality_measure_infty(x_k, g);
        g_norm = g0_norm;

        // print out - just to have idea how we are starting
        if(this->verbose_)
        {
          this->init_solver("QUasiTrustRegionVariableBound",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "ared","pred", "rho", "delta_k", "|| p_k ||_{inf} "});
          PrintInfo::print_iter_status(it, {g_norm});
        }

        it++;
        fun.value(x_k, E_old);


        UTOPIA_NO_ALLOC_BEGIN("Quasi TR Bound 1");
        // solve starts here
        while(!converged)
        {

          // ----------------------------------------------------------------------------
          //     new step p_k w.r. ||p_k|| <= delta
          // ----------------------------------------------------------------------------
          if(MatrixFreeQPSolver * tr_subproblem = dynamic_cast<MatrixFreeQPSolver*>(this->linear_solver().get()))
          {
            p_k.set(0.0);
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta);
            tr_subproblem->set_box_constraints(box);
            g_help = -1.0*g;
            tr_subproblem->solve(*multiplication_action, g_help, p_k);
            this->solution_status_.num_linear_solves++;
          }
          else
          {
            utopia_warning("QUasiTrustRegionVariableBound::Set suitable TR subproblem.... \n ");
          }

          pred = this->get_pred(g, *multiplication_action, p_k);
          // ----------------------------------------------------------------------------
          // ----------------------------------------------------------------------------
          // trial point
          x_trial = x_k + p_k;


          // value of the objective function with correction
          fun.value(x_trial, E_new);

          // decrease ratio
          ared = E_old - E_new;           // reduction observed on objective function
          pred = std::abs(pred);
          rho = ared/ pred;               // decrease ratio

    //----------------------------------------------------------------------------
    //     acceptance of trial point
    //----------------------------------------------------------------------------
          if(ared < 0.0 || pred < 0.0)
          {
            rho = 0.0;
          }
          else if(ared == pred)
          {
            rho = 1.0;
          }

          if (rho >= this->rho_tol()){
            it_successful_++;
          }


          // good reduction, accept trial point
          if (rho >= this->rho_tol())
          {
            x_k   = x_trial;
            E_old = E_new;

            y = g;
            fun.gradient(x_k, g);
            y = g - y;
          }
          // otherwise, keep old point
          else
          {
            fun.gradient(x_trial, g_help);
            y = g_help - g;
          }

          this->update(p_k, y, x_k, g);
    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          g_norm = this->criticality_measure_infty(x_k, g);
          s_norm = norm_infty(p_k);

          if(this->verbose_){
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, ared, pred, rho, delta, s_norm});
          }

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta);
    //----------------------------------------------------------------------------
    //      tr. radius update
    //----------------------------------------------------------------------------
          this->delta_update_inf(rho, p_k, delta);
          // this->delta_update_new(rho, p_k, delta, g);

          it++;
        }

          UTOPIA_NO_ALLOC_END();
          return false;
      }


    void set_trust_region_strategy(const std::shared_ptr<MatrixFreeQPSolver> &tr_linear_solver)
    {
      NonLinearSolver::set_linear_solver(tr_linear_solver);
    }


    private:
        void init_memory(const Layout &layout)
        {
            p_k.zeros(layout);
            g.zeros(layout);
            g_help.zeros(layout);
            y.zeros(layout);
            x_trial.zeros(layout);

            TrustRegionBase::init_memory(layout);

            initialized_ = true;
            layout_ = layout;
        }

        // void delta_update_new(const Scalar & rho, const Vector & p_k, Scalar & delta, const Vector & Pg)
        // {
        //   Vector Bp_k;

        //   auto multiplication_action = this->hessian_approx_strategy_->build_apply_H();
        //   multiplication_action->apply(p_k, Bp_k);

        //   Scalar contraction_factor = dot(Pg, p_k)/dot(p_k, Bp_k) * norm_infty(p_k);
        //   // std::cout<<"contraction_factor: "<< contraction_factor << "  \n";

        //   contraction_factor = std::abs(contraction_factor);

        //   Scalar b0 = 10e-5;
        //   Scalar b1 = 0.2;
        //   Scalar b2 = 0.8;

        //   Scalar sigma0 = 0.25;
        //   Scalar sigma1 = 0.5;
        //   Scalar sigma2 = 2.0;

        //   // std::cout<<"tho: "<< rho << "  \n";

        //   if(rho < b0)
        //   {
        //     delta = sigma0 * delta;
        //   }
        //   else if( b0 <= rho <= b1){
        //     delta = sigma1 * contraction_factor;
        //   }
        //   else if(b1 <= rho <= b2){
        //     delta = contraction_factor;
        //   }
        //   else{
        //     delta = sigma2 * contraction_factor;
        //   }

        //   // std::cout<<"delta: "<< delta << "  \n";


        // }


      SizeType it_successful_;
      Vector g, y, p_k, x_trial, g_help;
      bool initialized_;
      Layout layout_;

  };

}

#endif //UTOPIA_QUASI_TRUST_REGION_VARIABLE_CONSTRAINT_TR_HPP

