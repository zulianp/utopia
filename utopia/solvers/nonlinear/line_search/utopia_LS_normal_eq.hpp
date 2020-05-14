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
  class LeastSquaresNewton final: public NonLinearLeastSquaresSolver<Matrix, Vector>
  {
    typedef typename NonLinearLeastSquaresSolver<Matrix, Vector>::Solver Solver;
    using LSStrategy = utopia::LSStrategy<Vector>;
    using Scalar = typename utopia::Traits<Vector>::Scalar;
    using SizeType = typename utopia::Traits<Vector>::SizeType;

public:

    LeastSquaresNewton( const std::shared_ptr<Solver> &linear_solver= std::make_shared<ConjugateGradient<Matrix, Vector> >()):
                        NonLinearLeastSquaresSolver<Matrix, Vector>(linear_solver)
    {

    }

    bool solve(LeastSquaresFunction<Matrix, Vector> &fun, Vector &x_k)
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
        norms2(r_k, p_k, r_norm, s_norm);
        rel_norm = r_norm/r0_norm;

        // print iteration status on every iteration
        if(this->verbose_)
          PrintInfo::print_iter_status(it, {r_norm, rel_norm, s_norm, alpha_k});

              // check convergence and print interation info
        converged = this->check_convergence(it, r_norm, rel_norm, s_norm);

        it++;
      }
      return true;
    }


    void set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy)
    {
      ls_strategy_ = strategy;
    }

    void read(Input &in) override
    {
        NonLinearLeastSquaresSolver<Matrix, Vector>::read(in);
        if(ls_strategy_) {
            in.get("line-search", *ls_strategy_);
        }
    }


    void print_usage(std::ostream &os) const override
    {
        NonLinearLeastSquaresSolver<Matrix, Vector>::print_usage(os);
        this->print_param_usage(os, "line-search", "LSStrategy", "Input parameters for line-search strategy.", "-");
    }

  private:
    std::shared_ptr<LSStrategy> ls_strategy_;     /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */
  };

}

#endif //UTOPIA_LEAST_SQUARES_NEWTON_HPP

