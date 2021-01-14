#ifndef UTOPIA_SOLVER_NEWTON_HPP
#define UTOPIA_SOLVER_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>

namespace utopia {
/**
 * @brief      The Newton solver.
 *             Solver doesn't contain any globalization strategy, but it is
 * possible to set-up damping parameter.
 *
 *             The new iterate is obtained as:
 *              \f$ x_{k+1} = x_{k} - \alpha * H(x_k)^{-1} g(x_k), \f$ \n
 *              where \f$ \alpha  \f$ is duming parameter,
 *                    \f$ g  \f$ is gradient and \f$ H  \f$ is Hessian.
 *              Default value of \f$ \alpha  \f$ is 1, therefore full Newton
 * step.
 *
 * Example usage:
 * @snippet tests/utopia_SolverTest.cpp Newton CG example
 *
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
class Newton final : public NewtonBase<Matrix, Vector> {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
  using LSStrategy = utopia::LSStrategy<Vector>;

 public:
  Newton(const std::shared_ptr<Solver> &linear_solver =
             std::make_shared<ConjugateGradient<Matrix, Vector> >())
      : NewtonBase<Matrix, Vector>(linear_solver), alpha_(1.0) {}

  bool solve(Function<Matrix, Vector> &fun, Vector &x) override {
    using namespace utopia;

    if (empty(x)) {
      utopia_error("utopia::Newton, initial guess is empty vector");
    }

    init_memory(layout(x));

    Matrix hessian, preconditioner;

    Scalar g_norm = 1, g0_norm = 1, r_norm = 1, s_norm = 1;
    SizeType it = 0;

    bool converged = false;

    // notify listener
    fun.update(x);

    fun.gradient(x, grad_);
    g0_norm = norm2(grad_);
    g_norm = g0_norm;

    this->init_solver(
        "NEWTON", {" it. ", "|| g ||", "J", "r_norm", "|| p_k || ", "alpha_k"});

    if (this->verbose_) PrintInfo::print_iter_status(it, {g_norm, 1, 0});
    it++;

    while (!converged) {
      // find direction step
      step_.set(0.0);

      // setting up adaptive stopping criterium for linear solver
      if (this->has_forcing_strategy()) {
        if (auto *iterative_solver =
                dynamic_cast<IterativeSolver<Matrix, Vector> *>(
                    this->linear_solver_.get())) {
          iterative_solver->atol(this->estimate_ls_atol(g_norm, it));
        } else {
          utopia_error(
              "utopia::Newton::you can not use inexact Newton with exact "
              "linear solver. ");
        }
      }

      if (this->has_preconditioned_solver() && fun.has_preconditioner()) {
        fun.hessian(x, hessian, preconditioner);
        grad_neg_ = -1.0 * grad_;

        if (!this->check_values(it, fun, x, grad_, hessian)) return false;

        this->linear_solve(hessian, preconditioner, grad_neg_, step_);
      } else {
        fun.hessian(x, hessian);
        grad_neg_ = -1.0 * grad_;

        if (!this->check_values(it, fun, x, grad_, hessian)) return false;

        this->linear_solve(hessian, grad_neg_, step_);
      }

      if (ls_strategy_) {
        ls_strategy_->get_alpha(fun, grad_, x, step_, alpha_);
        x += alpha_ * step_;
      } else {
        // update x
        if (fabs(alpha_ - 1) < std::numeric_limits<Scalar>::epsilon()) {
          x += step_;
        } else {
          x += alpha_ * step_;
        }
      }

      // notify listener
      fun.update(x);
      fun.gradient(x, grad_);

      // norms needed for convergence check
      norms2(grad_, step_, g_norm, s_norm);
      r_norm = g_norm / g0_norm;

      Scalar J;
      fun.value(x, J);

      // // print iteration status on every iteration
      if (this->verbose_)
        PrintInfo::print_iter_status(it, {g_norm, J, r_norm, s_norm, alpha_});

      // // check convergence and print interation info
      converged = this->check_convergence(it, g_norm, r_norm, s_norm);
      it++;
    }

    this->print_statistics(it);

    return true;
  }

  void read(Input &in) override {
    NewtonBase<Matrix, Vector>::read(in);
    in.get("dumping", alpha_);

    if (ls_strategy_) {
      in.get("line-search", *ls_strategy_);
    }
  }

  void print_usage(std::ostream &os) const override {
    NewtonBase<Matrix, Vector>::print_usage(os);

    this->print_param_usage(os, "dumping", "real", "Default step size.", "1.0");
    this->print_param_usage(os, "line-search", "LSStrategy",
                            "Input parameters for line-search strategy.", "-");
  }

  void init_memory(const Layout &layout) {
    // init of linear solver
    NewtonBase<Matrix, Vector>::init_memory(layout);

    // init of vectors
    grad_neg_.zeros(layout);
    step_.zeros(layout);
    grad_.zeros(layout);

    // init of vectors
    if (ls_strategy_) ls_strategy_->init_memory(layout);
  }

  /**
   * @brief      Sets strategy for computing step-size.
   *
   * @param[in]  strategy  The line-search strategy.
   *
   * @return
   */
  void set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy) {
    ls_strategy_ = strategy;
  }

 private:
  Scalar alpha_;                            /*!< Dumping parameter. */
  std::shared_ptr<LSStrategy> ls_strategy_; /*!< Strategy used in order to
                                               obtain step \f$ \alpha_k \f$ */
  Vector grad_neg_, step_, grad_;
};

}  // namespace utopia
#endif  // UTOPIA_SOLVER_NEWTON_HPP
