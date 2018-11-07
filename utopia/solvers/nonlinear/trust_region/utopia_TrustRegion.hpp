#ifndef UTOPIA_SOLVER_TRUSTREGION_HPP
#define UTOPIA_SOLVER_TRUSTREGION_HPP

#include "utopia_NewtonBase.hpp"
#include "utopia_TRBase.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_Parameters.hpp"
#include "utopia_SteihaugToint.hpp"


 namespace utopia
 {
    	template<class Matrix, class Vector>
      /**
       * @brief      Base class for all TR solvers. Contains all general routines related to TR solvers.
       *             Design of class allows to provide different TR strategies in order to solve TR subproblem.
       */
     	class TrustRegion : public NewtonBase<Matrix, Vector>,
                          public TrustRegionBase<Matrix, Vector>
      {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem;

        typedef utopia::TrustRegionBase<Matrix, Vector> TrustRegionBase;
        typedef utopia::NewtonBase<Matrix, Vector> NonLinearSolver;

     	public:
      TrustRegion(const std::shared_ptr<TRSubproblem> &tr_subproblem = std::make_shared<SteihaugToint<Matrix, Vector>>(),
                  const Parameters params = Parameters())
                  : NonLinearSolver(tr_subproblem, params), it_successful_(0)
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

      /*!
      \details
                Determine, wheater problem is linear, or we need nonlinear solve
      @note
      \param fun          - function with evaluation routines
      \param H            - hessian
      \param g            - gradient
      \param p_N          - Newton step
      \param x_k          - current iterate
        */
      virtual bool linear_solution_check(

        Function<Matrix, Vector> &fun,
        Vector & g,
        const Matrix & H,
        Vector & p_N,
        Vector & x_k)
      {
        this->linear_solve(H, -1 * g, p_N);
        fun.gradient(x_k + p_N, g);
        Scalar g_norm = norm2(g);

        if(g_norm < 1e-7)
        {
          x_k += p_N;
          std::cout<<"To solve linear problem, TR solver is not really needed ...";
          return true;
        }

        x_k += p_N;
        return false;
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
         bool converged = false;
         NumericalTollerance<Scalar> tol(this->atol(), this->rtol(), this->stol());

         Scalar delta, product, ared, pred, rho, E, E_k, E_k1;

         SizeType it = 0;
         it_successful_ = 0;
         Scalar g_norm, g0_norm, r_norm, s_norm = std::numeric_limits<Scalar>::infinity();

         bool rad_flg = false;

         Vector g, p_CP = x_k, p_N = x_k, p_k = x_k, x_k1 = x_k;
         Matrix H;

         fun.hessian(x_k, H);
         fun.gradient(x_k, g);

          #define DEBUG_mode
          //  #define LS_check

        // TR delta initialization
        delta =  this->delta_init(x_k , this->delta0(), rad_flg);


        g0_norm = norm2(g);
        g_norm = g0_norm;

        // print out - just to have idea how we are starting
        #ifdef DEBUG_mode
          if(this->verbose_)
          {
              this->init_solver("TRUST_REGION_BASE",
                                {" it. ", "|| g ||", "r_norm", "<g, dx>", "J_k", "J_{k+1/2}", "J_{k+1}", "ared", "pred",
                                 "rho", "delta_k", "|| p_k || "});
            PrintInfo::print_iter_status(it, {g_norm});
          }

        #else
          if(this->verbose_)
          {
            this->init_solver("TRUST_REGION_BASE", {" it. ", "|| g ||", "r_norm", "J_k", "J_{k+1}", "rho", "delta_k", "|| p_k ||"});
            PrintInfo::print_iter_status(it, {g_norm});
          }
        #endif


        // found out if there is a linear solution - or start with the newton step
        #ifdef LS_check
          if(linear_solution_check(fun, g, H, p_N, x_k))
          {
            if(this->verbose_) {
              TrustRegionBase::check_convergence(*this, tol, this->max_it(), 0, 0, 0, 0, 0);
            }
            return true;
          }
          x_k1 = x_k;
        #endif

        it++;
        fun.value(x_k, E);
        fun.gradient(x_k, g);

        // solve starts here
        while(!converged)
        {
          fun.value(x_k, E_k);
          fun.hessian(x_k, H);
    //----------------------------------------------------------------------------
    //     new step p_k w.r. ||p_k|| <= delta
    //----------------------------------------------------------------------------
          if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->linear_solver_.get()))
            tr_subproblem->tr_constrained_solve(H, g, p_k, delta);

          pred = this->get_pred(g, H, p_k);
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

          if(it == 1 && rad_flg)
          {
            delta = norm2(p_k);
            delta *= 0.2;
          }

          // value of the objective function with correction
          fun.value(x_k + p_k, E_k1);
          product = dot(g, p_k);            // just to do tests

          // decrease ratio
          ared = E_k - E_k1;                // reduction observed on objective function
          pred = std::abs(pred);
          rho = ared/ pred;               // decrease ratio

    //----------------------------------------------------------------------------
    //     acceptance of trial point
    //----------------------------------------------------------------------------
          if(ared < 0 || pred < 0)
          {
            rho = 0;
          }
          else if(ared == pred)
          {
            rho = 1;
          }

          if (rho >= this->rho_tol())
            it_successful_++;

          this->trial_point_acceptance(rho, E, E_k, E_k1, p_k, x_k, x_k1);
    //----------------------------------------------------------------------------
    //    convergence check
    //----------------------------------------------------------------------------
          x_k = x_k1;
          fun.gradient(x_k, g);
          g_norm = norm2(g);
          r_norm = g_norm/g0_norm;
          s_norm = norm2(p_k);

          #ifdef DEBUG_mode
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, product, E, E_k, E_k1, ared, pred, rho, delta, s_norm});
          #else
            if(this->verbose_)
              PrintInfo::print_iter_status(it, {g_norm, r_norm, E, E_k1, rho, delta, s_norm});
          #endif

            converged = TrustRegionBase::check_convergence(*this, tol, this->max_it(), it, g_norm, r_norm, s_norm, delta);
    //----------------------------------------------------------------------------
    //      tr. radius update
    //----------------------------------------------------------------------------
          this->delta_update(rho, p_k, delta);
          it++;
        }

        // some benchmarking
        this->print_statistics(it);

          return true;
      }

      virtual void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector> > &ls) override
      {
          auto linear_solver = this->linear_solver();
          if (dynamic_cast<TRSubproblem *>(linear_solver.get()) != nullptr)
          {
              TRSubproblem * tr_sub = dynamic_cast<TRSubproblem *>(linear_solver.get());
              tr_sub->set_linear_solver(ls);
          }
      }

      virtual void set_trust_region_strategy(const std::shared_ptr<TRSubproblem> &tr_linear_solver)
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

#endif //UTOPIA_SOLVER_TRUSTREGION_HPP

