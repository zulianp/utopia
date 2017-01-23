/*
* @Author: alenakopanicakova
* @Date:   2016-05-11
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-11-07
*/

#ifndef UTOPIA_SOLVER_TRUSTREGION_BASE_HPP
#define UTOPIA_SOLVER_TRUSTREGION_BASE_HPP
#include "utopia_NonLinearSolver.hpp"
#include "utopia_TRSubproblem.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_Parameters.hpp"    
#include "utopia_NumericalTollerance.hpp"

namespace utopia 
{
    	template<class Matrix, class Vector>
      /**
       * @brief      Base class for all TR solvers. Contains all general routines related to TR solvers.
       *             Design of class allows to provide different TR strategies in order to solve TR subproblem. 
       */ 
  class TrustRegionBase 
  {
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
    typedef utopia::LinearSolver<Matrix, Vector> Solver;
    typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem; 

  public:
    TrustRegionBase(const Parameters params = Parameters())
    {

      set_parameters(params);        
    }

      /* @brief      Sets the parameters.
      *
      * @param[in]  params  The parameters
      */
    virtual void set_parameters(const Parameters params) 
    {
      delta_max_  = params.delta_max();
      delta0_     = params.delta0();
      gamma1_     = params.gamma1();
      gamma2_     = params.gamma2();
      eta1_       = params.eta1();
      eta2_       = params.eta2();
      rho_tol_    = params.rho_tol();
      eps_        = params.eps();

    }

    Scalar delta_max()  const  { return delta_max_; } 
    Scalar delta0()     const  { return delta0_; } 
    Scalar gamma1()     const  { return gamma1_; } 
    Scalar gamma2()     const  { return gamma2_; } 
    Scalar eta1()       const  { return eta1_; } 
    Scalar eta2()       const  { return eta2_; } 
    Scalar rho_tol()    const  { return rho_tol_; } 
    Scalar eps()        const  { return eps_; } 

    void delta_max(const bool & delta_max_in ) { delta_max_ = delta_max_in; }; 
    void delta0(const bool & delta0_in ) { delta0_ = delta0_in; }; 
    void gamma1(const bool & gamma1_in ) { gamma1_ = gamma1_in; }; 
    void gamma2(const bool & gamma2_in ) { gamma2_ = gamma2_in; }; 
    void eta1(const bool & eta1_in ) { eta1_ = eta1_in; }; 
    void eta2(const bool & eta2_in ) { eta2_ = eta2_in; }; 
    void rho_tol(const bool & rho_tol_in ) { rho_tol_ = rho_tol_in; }; 
    void eps(const bool & eps_in ) { eps_ = eps_in; }; 

  protected:

    /**
     * @brief      Calculates the predicate reduction  m_k(0) - m_k(p_k)
     *
     * @param[in]  g     Gradient
     * @param[in]  H     Hessian
     * @param[in]  p_k   current step
     * @param      pred  predicted reduction
     */
    virtual void compute_pred_red( const Vector & g, const Matrix & H, const Vector & p_k, Scalar &pred)
    {
    	Scalar l_term = dot(g, p_k);
    	Scalar qp_term = dot(p_k, H * p_k);
    	pred = - l_term - 0.5 * qp_term; 

    }


    virtual bool check_convergence(
      Monitor<Matrix, Vector> &monitor,
      const NumericalTollerance<Scalar> &tol,
      const SizeType max_it,
      const SizeType &it, 
      const Scalar & g_norm, 
      const Scalar & r_norm, 
      const Scalar & s_norm, 
      const Scalar & delta) const
    {   
        // termination because norm of grad is down
      if(g_norm < tol.absolute_tollerance())
      {
        monitor.exit_solver(it, ConvergenceReason::CONVERGED_FNORM_ABS);
        return true; 
      }

        // step size so small that we rather exit than wait for nan's
      if(s_norm < tol.step_tollerance())
      {
        monitor.exit_solver(it, ConvergenceReason::CONVERGED_SNORM_RELATIVE);
        return true; 
      }

        // step size so small that we rather exit than wait for nan's
      if(r_norm < tol.relative_tollerance())
      {
        monitor.exit_solver(it, ConvergenceReason::CONVERGED_FNORM_RELATIVE);
        return true; 
      }

        // check number of iterations
      if( it > max_it)
      {
        monitor.exit_solver(it, ConvergenceReason::DIVERGED_MAX_IT);
        return true; 
      }

        // do not hard code this 
      if(delta < eps_)
      {
        monitor.exit_solver(it, ConvergenceReason::CONVERGED_TR_DELTA);
        return true; 
      }

      return false; 
    }


    /*!
    \details
              Trial point acceptance 
    @note
    \param rho            - ared/pred
    \param E              - Energy after acceptance
    \param E_k            - current energy
    \param E_k1           - next iterate energy
    \param p_k            - current step
    \param x_k            - curren iterate
    \param x_k1           - new iterate
      */
    virtual bool trial_point_acceptance(const Scalar &rho, Scalar &E, const Scalar &E_k, const Scalar &E_k1,  const Vector & p_k, const Vector & x_k, Vector & x_k1)
    {
      // good reduction, accept trial point 
      if (rho >= rho_tol_)
      {
        x_k1 = x_k + p_k;
        E = E_k1; 
      }
      // otherwise, keep old point
      else
      {
        x_k1 = x_k;
        E = E_k; 
      }

      return true; 
    }



    /*!
    \details
    TR radius update function
    @note
    \param rho            - ared/pred
    \param radius          - tr. radius
    \param p_k            - iterate step
      */
    virtual bool delta_update(const Scalar &rho, const Vector &p_k, Scalar &radius)
    {
      if(rho < eta1_)
      {
          // in case of L2 norm - this is good option 
          // for sloving with GS or MG in _inf norm, change ... 
        radius = gamma1_ * radius; 
      }

      else if (rho > eta2_ && (norm2(p_k) > ((1 - eps_) * radius ) ))
      {
        radius = std::min(gamma2_ * radius, delta_max_); 
      }      
      return true; 

    }


    /*!
    \details
    TR radius initialization
    - returns false - choice of tr radius was given by user, alg respect it and doesn't change
  
    @note
    \param x_k                - initial guess/ current iterate
    \param radius              - tr. radius 
      */
    virtual bool delta_init(const Vector & x_k , Scalar & radius)
    {
      Scalar x_norm = norm2(x_k); 
      if(radius == 0)
      { 
        radius =  (x_norm == 0)? 100 : 100 * x_norm; 
      }
      else
      {
        return false;  
      }

      return true; 
    }


    /**
     * @brief      Gets the prediction reduction for 
     *
     * @param[in]  g     gradient
     * @param[in]  B     Hessian
     * @param[in]  p_k    step 
     * @param      pred  The predicted reduction. 
     */
    virtual bool get_pred(const Vector & g, const Matrix & B, const Vector & p_k, Scalar &pred)
    {
      Scalar l_term = dot(g, p_k);
      Scalar qp_term = dot(p_k, B * p_k);
      pred = - l_term - 0.5 * qp_term; 
      return true; 
    }

  private: 
    Scalar delta_max_; 
    Scalar delta0_; 
    Scalar gamma1_;
    Scalar gamma2_;
    Scalar eta1_;
    Scalar eta2_;
    Scalar rho_tol_;
    Scalar eps_; 
  };

}

#endif //UTOPIA_SOLVER_TRUSTREGION_BASE_HPP

