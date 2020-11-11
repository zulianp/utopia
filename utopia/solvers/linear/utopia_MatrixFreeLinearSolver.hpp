#ifndef UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
#define UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP

#include "utopia_Instance.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Preconditioner.hpp"
#include "utopia_SolverForwardDeclarations.hpp"

namespace utopia {

template <class Vector>
class MatrixFreeLinearSolver : virtual public Configurable,
                               virtual public Clonable,
                               virtual public Preconditioner<Vector> {
 public:
  using Preconditioner<Vector>::init_memory;

  ~MatrixFreeLinearSolver() override = default;
  virtual bool solve(const Operator<Vector> &A, const Vector &rhs,
                     Vector &sol) = 0;

  /*! @brief if overriden the subclass has to also call this one first
   */
  void update(const Operator<Vector> & /*A*/) override =
      0;  //{ UTOPIA_UNUSED(A); }

  MatrixFreeLinearSolver *clone() const override = 0;

  void read(Input & /*in*/) override {}
  void print_usage(std::ostream & /*os*/) const override {}
};

template <class Matrix, class Vector>
class OperatorBasedLinearSolver : public MatrixFreeLinearSolver<Vector>,
                                  public PreconditionedSolver<Matrix, Vector> {
 public:
  using MatrixFreeLinearSolver<Vector>::update;
  using PreconditionedSolver<Matrix, Vector>::update;
  using MatrixFreeLinearSolver<Vector>::solve;

  ~OperatorBasedLinearSolver() override = default;

  bool solve(const Matrix &A, const Vector &b, Vector &x) override {
    update(make_ref(A));
    return solve(operator_cast<Vector>(A), b, x);
  }

  void update(const std::shared_ptr<const Matrix> &op) override {
    PreconditionedSolver<Matrix, Vector>::update(op);
    update(operator_cast<Vector>(*op));
  }

  bool smooth(const Vector &rhs, Vector &x) override {
    SizeType temp = this->max_it();
    this->max_it(this->sweeps());
    solve(operator_cast<Vector>(*this->get_operator()), rhs, x);
    this->max_it(temp);
    return true;
  }

  /**
   * @brief      Solution routine after update.
   *
   * @param[in]  b     The right hand side.
   * @param      x     The initial guess/solution.
   *
   * @return true if the linear system has been solved up to required
   * tollerance. False otherwise
   */
  bool apply(const Vector &b, Vector &x) override {
    return solve(operator_cast<Vector>(*this->get_operator()), b, x);
  }

  OperatorBasedLinearSolver *clone() const override = 0;

  void read(Input &in) override {
    MatrixFreeLinearSolver<Vector>::read(in);
    PreconditionedSolver<Matrix, Vector>::read(in);
  }

  void print_usage(std::ostream &os) const override {
    MatrixFreeLinearSolver<Vector>::print_usage(os);
    PreconditionedSolver<Matrix, Vector>::print_usage(os);
  }
};

template <class Vector>
class EmptyPrecondMatrixFreeLinearSolver final
    : public MatrixFreeLinearSolver<Vector> {
 public:
  void set_preconditioner(
      const std::shared_ptr<Preconditioner<Vector> > &precond) {
    precond_ = precond;
  }

  bool solve(const Operator<Vector> & /*A*/, const Vector &rhs,
             Vector &sol) override {
    return apply(rhs, sol);
  }

  bool apply(const Vector &rhs, Vector &sol) override {
    if (precond_) {
      precond_->apply(rhs, sol);
    } else {
      utopia_warning(
          "EmptyPrecondMatrixFreeLinearSolver: preconditioner is missing \n");
    }
    return true;
  }

  EmptyPrecondMatrixFreeLinearSolver *clone() const override {
    return new EmptyPrecondMatrixFreeLinearSolver(*this);
  }

  void read(Input &in) override {
    MatrixFreeLinearSolver<Vector>::read(in);
    if (precond_) {
      in.get("precond", *precond_);
    }
  }

  void update(const Operator<Vector> &A) override {
    if (precond_) {
      precond_->update(A);
    }
  }

  void print_usage(std::ostream &os) const override {
    MatrixFreeLinearSolver<Vector>::print_usage(os);
    this->print_param_usage(os, "precond", "Preconditioner",
                            "Input parameters for preconditioner.", "-");
  }

 private:
  std::shared_ptr<Preconditioner<Vector> > precond_;
};

}  // namespace utopia

#endif  // UTOPIA_MATRIX_FREE_LINEAR_SOLVER_HPP
