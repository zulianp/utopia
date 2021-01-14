#ifndef UTOPIA_PRECONDITIONED_SOLVER_INTERFACE_HPP
#define UTOPIA_PRECONDITIONED_SOLVER_INTERFACE_HPP

#include "utopia_IterativeSolver.hpp"
#include "utopia_LinearSolver.hpp"

namespace utopia {
template <class Vector>
class PreconditionedSolverInterface : virtual public Configurable,
                                      virtual public Clonable {
 public:
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  using Preconditioner = utopia::Preconditioner<Vector>;
  // typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
  // typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

  /**
   * @brief      Sets the preconditioner.
   *
   * @param[in]  precond  The precondition
   */
  virtual void set_preconditioner(
      const std::shared_ptr<Preconditioner> &precond) {
    precond_ = precond;
  }

  virtual std::shared_ptr<Preconditioner> &get_preconditioner() {
    return precond_;
  }

  /**
   * @brief      Resets the preconditioner.
   */
  virtual void reset_preconditioner() { precond_.reset(); }

  virtual const std::shared_ptr<Preconditioner> &get_preconditioner() const {
    return precond_;
  }

  virtual bool has_preconditioner() const {
    return static_cast<bool>(precond_);
  }

  virtual void update(const std::shared_ptr<Operator<Vector>> &op) {
    precond_->update(*op);
  }

  virtual void update(const std::shared_ptr<const Operator<Vector>> &op) {
    precond_->update(*op);
  }

  void read(Input &in) override {
    if (precond_) {
      in.get("precond", *precond_);
    }
  }

  void print_usage(std::ostream &os) const override {
    this->print_param_usage(os, "precond", "Preconditioner",
                            "Input parameters for preconditioner.", "-");
  }

  virtual PreconditionedSolverInterface &operator=(
      const PreconditionedSolverInterface &other) {
    if (this == &other) {
      return *this;
    }

    copy_preconditioner_from(other);
    return *this;
  }

  PreconditionedSolverInterface(const PreconditionedSolverInterface &other) {
    copy_preconditioner_from(other);
  }

  PreconditionedSolverInterface() {}

  virtual ~PreconditionedSolverInterface() {}

 protected:
  std::shared_ptr<Preconditioner> precond_;

  virtual void copy_preconditioner_from(
      const PreconditionedSolverInterface &other) {
    if (other.has_preconditioner()) {
      this->set_preconditioner(
          std::shared_ptr<Preconditioner>(other.precond_->clone()));
    }
  }
};
}  // namespace utopia

#endif  // UTOPIA_PRECONDITIONED_SOLVER_INTERFACE_HPP
