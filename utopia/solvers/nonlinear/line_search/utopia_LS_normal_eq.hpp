/*
* @Author: alenakopanicakova
* @Date:   2016-05-10
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-14
*/

#ifndef UTOPIA_LEAST_SQUARES_NEWTON_HPP
#define UTOPIA_LEAST_SQUARES_NEWTON_HPP

#include "utopia_LS_Strategy.hpp"
#include "utopia_SimpleBacktracking.hpp"     

namespace utopia 
{
      /**
       * @brief      Newton Line Search class specified for nonlinear least square eq. \n
       *
       */
  template<class Matrix, class Vector>
  class LeastSquaresNewton : public NonLinearLeastSquaresSolver<Matrix, Vector> 
  {
    typedef typename NonLinearLeastSquaresSolver<Matrix, Vector>::Solver Solver;
    typedef utopia::LSStrategy<Vector> LSStrategy; 
    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

  public:

    LeastSquaresNewton( const std::shared_ptr<Solver> &linear_solver= std::make_shared<ConjugateGradient<Matrix, Vector> >(), 
                        const Parameters params               = Parameters()): 
                        NonLinearLeastSquaresSolver<Matrix, Vector>(linear_solver, params) 
    {

    }

    virtual bool solve(LeastSquaresFunction<Matrix, Vector> &fun, Vector &x_k) 
    {

      if(!this->ls_strategy_)
        std::cerr <<"utopia:: Line-search base:: Missing ls_strategy. Please provide before executing solve. \n"; 

      using namespace utopia;
      this->init_solver("NEWTON LINE SEARCH - normal eq.",
        {" it. ", "|| r ||", "rel_norm", "|| p_k || ", "alpha_k"});

      Vector r_k, p_k, g; 
      Matrix J_k, J_T; 

      bool converged = false; 

      Scalar alpha_k, r_norm, r0_norm, rel_norm, s_norm = std::numeric_limits<Scalar>::max();
      SizeType it = 0; 

      fun.residual(x_k, r_k);
      r0_norm = norm2(r_k);
      r_norm = r0_norm;

      if(this->verbose_)
        PrintInfo::print_iter_status(it, {r_norm, 1, 0, 1}); 
      it++; 

      while(!converged)
      {  
        fun.jacobian(x_k, J_k); 
        J_T = transpose(J_k); 

              //  Gauss-Newton step
        this->linear_solve(J_k, -1 * r_k, p_k);

              // this sould go away in future 
        g = J_T * r_k; 

        if(this->ls_strategy_)
          this->ls_strategy_->get_alpha(fun, g, x_k, p_k, alpha_k); 
        
        x_k += alpha_k * p_k; 


        fun.residual(x_k, r_k);

              // norms needed for convergence check 
        r_norm = norm2(r_k);
        rel_norm = r_norm/r0_norm; 
        s_norm = norm2(p_k); 

              // print iteration status on every iteration 
        if(this->verbose_)
          PrintInfo::print_iter_status(it, {r_norm, rel_norm, s_norm, alpha_k}); 

              // check convergence and print interation info
        converged = this->check_convergence(it, r_norm, rel_norm, s_norm); 

        it++; 
      } 
      return true;
    }


    virtual bool set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy)
    {
      ls_strategy_ = strategy; 
      ls_strategy_->set_parameters(this->parameters());
      return true; 
    }

  private:
    std::shared_ptr<LSStrategy> ls_strategy_;     /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */  
  };

}

#endif //UTOPIA_LEAST_SQUARES_NEWTON_HPP

