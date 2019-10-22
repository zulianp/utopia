#ifndef UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP
#define UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP

#include "utopia_NewtonBase.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_Allocations.hpp"

 namespace utopia
 {
        template<class Matrix, class Vector>
         class TrustRegionVariableBound final: public VariableBoundSolverInterface<Vector>,
                                            public TrustRegionBase<Vector>,
                                            public NewtonBase<Matrix, Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::QPSolver<Matrix, Vector>      QPSolver;

        typedef utopia::TrustRegionBase<Vector>       TrustRegionBase;
        typedef utopia::NewtonBase<Matrix, Vector>    NewtonBase;

         public:
      TrustRegionVariableBound( const std::shared_ptr<QPSolver> &tr_subproblem) :
                                NewtonBase(tr_subproblem)
      {

      }

      using TrustRegionBase::get_pred;

      void read(Input &in) override
      {
        TrustRegionBase::read(in);
        NewtonBase::read(in);
      }

      void print_usage(std::ostream &os) const override
      {
        TrustRegionBase::print_usage(os);
        NewtonBase::print_usage(os);
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
        SizeType it_successful = 0;
        const Scalar infty = std::numeric_limits<Scalar>::infinity();
        Scalar g_norm = infty, g0_norm = infty, r_norm = infty, s_norm = infty;
        bool rad_flg = false;

        Vector g = 0*x_k, p_k = 0*x_k, x_k1 = 0*x_k;
        Matrix H;

        this->make_iterate_feasible(x_k);
        fun.gradient(x_k, g);

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);
        // delta = 10;        // testing

        g0_norm = this->criticality_measure_infty(x_k, g);
        g_norm = g0_norm;

        // print out - just to have idea how we are starting
        if(this->verbose_)
        {
          this->init_solver("TRUST_REGION_BASE",
                              {" it. ", "||P_c(x-g)-x||","J_k", "J_{k+1}", "ared","pred", "rho", "delta_k", "|| p_k || "});
          PrintInfo::print_iter_status(it, {g_norm});
        }

        it++;
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
          if(QPSolver * tr_subproblem = dynamic_cast<QPSolver*>(this->linear_solver_.get()))
          {
            // UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound0");
            p_k.set(0.0); 
            auto box = this->merge_pointwise_constraints_with_uniform_bounds(x_k, -1.0 * delta, delta);
            // UTOPIA_NO_ALLOC_END();
            
            UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound1");
            tr_subproblem->set_box_constraints(box);
            
            if(accepted){
              g = -g; 
            }
            
            UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound2");
            tr_subproblem->solve(H, g, p_k);
            // UTOPIA_NO_ALLOC_END();
            this->solution_status_.num_linear_solves++;
          }
          else
          {
            utopia_warning("TrustRegionVariableBound::Set suitable TR subproblem.... \n ");
          }

          UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound3");
          pred = this->get_pred(g, H, p_k);
          UTOPIA_NO_ALLOC_END();
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

          UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound5");
          accepted = this->trial_point_acceptance(rho, x_k1, x_k);
          UTOPIA_NO_ALLOC_END();

          if (rho >= this->rho_tol())
            it_successful++;

    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          if(accepted)
          {
            fun.gradient(x_k, g);
            UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound6");
            g_norm = this->criticality_measure_infty(x_k, g);
            UTOPIA_NO_ALLOC_END();
          }

          UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound7");
          s_norm = norm_infty(p_k);
          UTOPIA_NO_ALLOC_END();

          if(this->verbose_)
            PrintInfo::print_iter_status(it, {g_norm, E_old, E_new, ared, pred, rho, delta, s_norm});

          converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta);

    //----------------------------------------------------------------------------
    //      tr. radius update
    //----------------------------------------------------------------------------
          UTOPIA_NO_ALLOC_BEGIN("TrustRegionVariableBound8");
          this->delta_update(rho, p_k, delta, true);
          UTOPIA_NO_ALLOC_END();
          it++;
        }

        // some benchmarking
        TrustRegionBase::print_statistics(it, it_successful);

          return false;
      }


      void set_trust_region_strategy(const std::shared_ptr<QPSolver> &tr_linear_solver)
      {
        NewtonBase::set_linear_solver(tr_linear_solver);
      }

    private:
      Scalar get_pred(const Vector & g, const Matrix & B, const Vector & p_k)
      {
        B_pk_ = B * p_k; 
        return (dot(g, p_k) -0.5 *dot(B_pk_, p_k));
      }


      Vector B_pk_; 

  };

}

#endif //UTOPIA_SOLVER_BOX_CONSTRAINT_TR_HPP

