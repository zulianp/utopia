#ifndef UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
#define UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP

#include <utility>

#include "utopia_ConjugateGradient.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_InexactNewtonInterface.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_PrintInfo.hpp"

namespace utopia {

template <class Matrix, class Vector>
class NewtonBase : public NonLinearSolver<Vector>,
                   public InexactNewtonInterface<Vector> {
 public:
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  using Solver = utopia::LinearSolver<Matrix, Vector>;
  using DiffController = utopia::DiffController<Matrix, Vector>;

  NewtonBase(std::shared_ptr<Solver> linear_solver)
      : NonLinearSolver<Vector>(),
        InexactNewtonInterface<Vector>(),
        linear_solver_(std::move(linear_solver)),
        check_diff_(false) {}

  ~NewtonBase() override = default;

  virtual bool solve(Function<Matrix, Vector> &fun, Vector &x) = 0;

  virtual bool solve(Function_rhs<Matrix, Vector> &fun, Vector &x,
                     const Vector &rhs) {
    fun.set_rhs(rhs);
    bool converged = this->solve(fun, x);
    fun.reset_rhs();
    return converged;
  }

  /**
   * @brief      Enables the differentiation control.
   *
   * @param[in]  checkDiff  Option, if eanable diff_control or no.
   */
  void enable_differentiation_control(bool checkDiff) {
    check_diff_ = checkDiff;
  }

  inline bool differentiation_control_enabled() const { return check_diff_; }

  bool check_values(const SizeType iterations,
                    const Function<Matrix, Vector> &fun, const Vector &x,
                    const Vector &gradient, const Matrix &hessian) {
    if (check_diff_ && !controller_.check(fun, x, gradient, hessian)) {
      this->exit_solver(iterations, ConvergenceReason::DIVERGED_INNER);
      return false;
    }

    return true;
  }

  void read(Input &in) override {
    NonLinearSolver<Vector>::read(in);
    in.get("check-diff", check_diff_);

    if (check_diff_) {
      in.get("diff-controller", controller_);
    }

    if (linear_solver_) {
      in.get("linear-solver", *linear_solver_);
    }
  }

  void print_usage(std::ostream &os) const override {
    NonLinearSolver<Vector>::print_usage(os);
    this->print_param_usage(os, "check_diff", "bool",
                            "Enables finite difference controller",
                            std::to_string(check_diff_));

    if (linear_solver_) {
      this->print_param_usage(os, "linear-solver", "LinearSolver",
                              "Input parameters for linear solver.", "-");
      linear_solver_->print_usage(os);
    } else {
      this->print_param_usage(os, "linear-solver", "LinearSolver",
                              "Input parameters for linear solver.",
                              "- (null)");
    }
  }

  /**
   * @brief      Changes linear solver used inside of nonlinear-solver.
   *
   * @param[in]  linear_solver  The linear solver
   */
  virtual void set_linear_solver(const std::shared_ptr<Solver> &linear_solver) {
    linear_solver_ = linear_solver;
  }

  inline DiffController &controller() { return controller_; }

 public:
  inline std::shared_ptr<Solver> linear_solver() const {
    return linear_solver_;
  }

 protected:
  inline bool linear_solve(const Matrix &mat, const Vector &rhs, Vector &sol) {
    linear_solver_->update(make_ref(mat));
    this->solution_status_.num_linear_solves++;
    return linear_solver_->apply(rhs, sol);
  }

  inline bool has_preconditioned_solver() {
    return dynamic_cast<PreconditionedSolver<Matrix, Vector> *>(
        linear_solver_.get());
  }

  inline bool linear_solve(const Matrix &mat, const Matrix &prec,
                           const Vector &rhs, Vector &sol) {
    static_cast<PreconditionedSolver<Matrix, Vector> *>(linear_solver_.get())
        ->update(make_ref(mat), make_ref(prec));
    this->solution_status_.num_linear_solves++;
    return linear_solver_->apply(rhs, sol);
  }

  void init_memory(const Layout &layout) {
    linear_solver_->init_memory(layout);
  }

  std::shared_ptr<Solver> linear_solver_; /*!< Linear solver parameters. */
  DiffController controller_;
  bool check_diff_; /*!< Enable differentiation control. */
};
}  // namespace utopia

#endif  // UTOPIA_NEWTON_BASED_NONLINEAR_SOLVER_HPP
