#ifndef UTOPIA_CHEBYSHEV_HPP
#define UTOPIA_CHEBYSHEV_HPP

#include <memory>
#include "utopia_MatrixFreeLinearSolver.hpp"

namespace utopia {

/**
 * @brief      Conjugate Gradient solver. Works with all utopia tensor types.
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
class Chebyshev3level final : public OperatorBasedLinearSolver<Matrix, Vector> {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  typedef utopia::LinearSolver<Matrix, Vector> Solver;
  using Preconditioner = utopia::Preconditioner<Vector>;

 public:
  using Super = utopia::OperatorBasedLinearSolver<Matrix, Vector>;
  using Super::solve;
  using Super::update;

  Chebyshev3level() : eps_eig_est_(1e-5), power_method_max_it_(5) {}

  Chebyshev3level(const Chebyshev3level &other)
      : eps_eig_est_(1e-5), power_method_max_it_(5) {}

  void read(Input &in) override {
    OperatorBasedLinearSolver<Matrix, Vector>::read(in);
    in.get("eig_comp_tol", eps_eig_est_);
    in.get("power_method_max_it", power_method_max_it_);
  }

  void init_memory(const Layout &layout) override {
    assert(layout.local_size() > 0);
    OperatorBasedLinearSolver<Matrix, Vector>::init_memory(layout);

    // resets all buffers in case the size has changed
    r_.zeros(layout);
    help_f1.zeros(layout);
    help_f2.zeros(layout);
    fi_.zeros(layout);
    p_.zeros(layout);

    initialized_ = true;
    layout_ = layout;
  }

  void print_usage(std::ostream &os) const override {
    OperatorBasedLinearSolver<Matrix, Vector>::print_usage(os);
  }

  void update(const Operator<Vector> &A) override {
    const auto layout_rhs = row_layout(A);

    if (!initialized_ || !layout_rhs.same(layout_)) {
      init_memory(layout_rhs);
    }
  }

  Chebyshev3level *clone() const override { return new Chebyshev3level(*this); }

  void copy(const Chebyshev3level &other) {
    Super::operator=(other);
    // reset_initial_guess_ = other.reset_initial_guess_;
  }

  bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) override {
    // std::cout << "----------- here ------------- \n";
    Scalar eigMax = 1.1 * this->get_max_eig(A);
    Scalar eigMin = eigMax / 30;

    Scalar avg_eig = (eigMax + eigMin) / 2.0;
    Scalar diff_eig = (eigMax - eigMin) / 2.0;

    A.apply(x, help_f2);
    r_ = help_f2 - b;

    SizeType it = 0;
    Scalar r_norm = 9e9, alpha = 0.0, beta = 0.0;
    bool converged = false;

    while (!converged) {
      if (it == 0) {
        p_ = -1.0 * r_;
        alpha = 1.0 / avg_eig;
      } else if (it == 1) {
        beta = 0.5 * (diff_eig * alpha) * (diff_eig * alpha);
        alpha = 1.0 / (avg_eig - beta);
        p_ = -1.0 * r_ + (beta * p_);
      } else {
        beta = 0.25 * (diff_eig * alpha) * (diff_eig * alpha);
        alpha = 1.0 / (avg_eig - beta);
        p_ = -1.0 * r_ + (beta * p_);
      }

      x = x + (alpha * p_);
      A.apply(p_, help_f2);
      r_ = r_ + (alpha * help_f2);

      r_norm = norm2(r_);

      if (this->verbose()) {
        PrintInfo::print_iter_status(it, {r_norm});
      }

      converged = this->check_convergence(it, r_norm, 1, 1);
      it++;
    }

    return true;
  }

 private:
  Scalar get_max_eig(const Operator<Vector> &A) {
    // Super simple power method to estimate the biggest eigenvalue
    assert(!empty(help_f2));

    // TODO:: setup random vector
    // help_f2.set(1.0);
    {
      auto d_help_f2 = local_view_device(help_f2);

      parallel_for(local_range_device(help_f2),
                   UTOPIA_LAMBDA(const SizeType i) {

                     const Scalar val = ((Scalar)std::rand() / (RAND_MAX)) + 1;

                     d_help_f2.set(i, val);
                   });
    }

    SizeType it = 0;
    bool converged = false;
    Scalar gnorm, lambda = 0.0, lambda_old;

    while (!converged) {
      help_f1 = help_f2;
      A.apply(help_f1, help_f2);
      help_f2 = Scalar(1.0 / Scalar(norm2(help_f2))) * help_f2;

      lambda_old = lambda;

      A.apply(help_f2, help_f1);
      lambda = dot(help_f2, help_f1);

      fi_ = help_f2 - help_f1;
      gnorm = norm2(fi_);

      converged = ((gnorm < eps_eig_est_) ||
                   (std::abs(lambda_old - lambda) < eps_eig_est_) ||
                   it > power_method_max_it_)
                      ? true
                      : false;

      it = it + 1;
    }

    if (this->verbose())
      utopia::out() << "Power method converged in " << it
                    << " iterations. Largest eig: " << lambda << "  \n";

    return lambda;
  }

  bool initialized_{false};
  Layout layout_;
  Scalar eps_eig_est_;
  SizeType power_method_max_it_;

  // This fields are not to be copied anywhere
  Vector r_, p_, help_f1, help_f2, fi_;
};
}  // namespace utopia

#endif
