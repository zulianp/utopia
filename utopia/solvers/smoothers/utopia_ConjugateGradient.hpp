#ifndef UTOPIA_CONJUGATE_GRADIENT_HPP
#define UTOPIA_CONJUGATE_GRADIENT_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

/**
 * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
class ConjugateGradient final
    : public OperatorBasedLinearSolver<Matrix, Vector> {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  typedef utopia::LinearSolver<Matrix, Vector> Solver;
  using Preconditioner = utopia::Preconditioner<Vector>;

 public:
  using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
  using Super::solve;
  using Super::update;

  ConjugateGradient();

  void reset_initial_guess(const bool val);
  inline void apply_gradient_descent_step(const bool val) {
    apply_gradient_descent_step_ = val;
  }

  void read(Input &in) override;

  void init_memory(const Layout &layout) override;

  void print_usage(std::ostream &os) const override;

  bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override;

  void update(const Operator<Vector> &A) override;

  ConjugateGradient *clone() const override;

  void copy(const ConjugateGradient &other);
  ConjugateGradient(const ConjugateGradient &other);

 private:
  bool unpreconditioned_solve(const Operator<Vector> &A, const Vector &b,
                              Vector &x);

  bool preconditioned_solve(const Operator<Vector> &A, const Vector &b,
                            Vector &x);

  bool check_solution(const Operator<Vector> &A, const Vector &x,
                      const Vector &b) const;

  void gradient_descent_step(const Operator<Vector> &A, const Vector &b,
                             Vector &x);

  bool reset_initial_guess_{false};
  bool initialized_{false};
  bool apply_gradient_descent_step_{false};
  Layout layout_;

  // This fields are not to be copied anywhere
  Vector r, p, q, Ap, r_new, z, z_new;
};
}  // namespace utopia

#endif  // UTOPIA_CONJUGATE_GRADIENT_HPP
