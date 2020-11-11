#ifndef UTOPIA_SOLVER_QUASI_NEWTON_HPP
#define UTOPIA_SOLVER_QUASI_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_Layout.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_QuasiNewtonBase.hpp"

#include <iomanip>
#include <limits>

namespace utopia {

template <class Matrix, class Vector>
class QuasiNewton : public QuasiNewtonBase<Vector> {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  using HessianApproximation = utopia::HessianApproximation<Vector>;
  using LinSolver = utopia::MatrixFreeLinearSolver<Vector>;

 public:
  using utopia::QuasiNewtonBase<Vector>::init_memory;

  QuasiNewton(const std::shared_ptr<HessianApproximation> &hessian_approx,
              const std::shared_ptr<LinSolver> &linear_solver)
      : QuasiNewtonBase<Vector>(hessian_approx, linear_solver),
        initialized_(false) {}

  QuasiNewton(const std::shared_ptr<LinSolver> &linear_solver)
      : QuasiNewtonBase<Vector>(linear_solver), initialized_(false) {}

  bool solve(FunctionBase<Vector> &fun, Vector &x) override {
    using namespace utopia;

    Scalar g_norm, g0_norm, r_norm = 1, s_norm = 1;
    SizeType it = 0;
    Scalar alpha = 1.0;

    bool converged = false;

    auto x_layout = utopia::layout(x);
    if (!initialized_ || !x_layout.same(layout_)) {
      init_memory(x_layout);
    }

    fun.gradient(x, g);
    g0_norm = norm2(g);
    g_norm = g0_norm;

    QuasiNewtonBase<Vector>::init_memory(x, g);

    if (this->verbose_) {
      this->init_solver("QUASI NEWTON",
                        {" it. ", "|| g ||", "r_norm", "|| p_k || ", "alpha"});
      PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm});
    }
    it++;

    UTOPIA_NO_ALLOC_BEGIN("Quasi_Newton");
    while (!converged) {
      // UTOPIA_NO_ALLOC_BEGIN("Quasi1");
      g_minus = -1.0 * g;
      this->linear_solve(g_minus, s);
      // UTOPIA_NO_ALLOC_END();

      // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:2");
      alpha = this->get_alpha(fun, g, x, s);
      // UTOPIA_NO_ALLOC_END();

      // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:2.1");
      s *= alpha;
      x += s;
      // UTOPIA_NO_ALLOC_END();

      // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:3");
      y = g;
      // UTOPIA_NO_ALLOC_END();
      fun.gradient(x, g);

      // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:3.1");
      // norms needed for convergence check
      norms2(g, s, g_norm, s_norm);
      r_norm = g_norm / g0_norm;

      // diff between fresh and old grad...
      y = g - y;
      // UTOPIA_NO_ALLOC_END();

      // UTOPIA_NO_ALLOC_BEGIN("Quasi Newton:4");
      this->update(s, y, x, g);
      // UTOPIA_NO_ALLOC_END();

      // print iteration status on every iteration
      if (this->verbose_)
        PrintInfo::print_iter_status(it, {g_norm, r_norm, s_norm, alpha});

      // check convergence and print interation info
      converged = this->check_convergence(it, g_norm, r_norm, s_norm);

      it++;
    }
    UTOPIA_NO_ALLOC_END();

    this->print_statistics(it);
    return true;
  }

  void update(const Vector &s, const Vector &y, const Vector &x,
              const Vector &g) override {
    std::cout << "------ change this ------- \n";
    if (auto *JFNK_mg = dynamic_cast<JFNK_Multigrid<Matrix, Vector> *>(
            this->mf_linear_solver_.get())) {
      JFNK_mg->update(s, y, x, g);
    } else {
      QuasiNewtonBase<Vector>::update(s, y, x, g);
    }
  }

  void initialize_approximation(const Vector &x, const Vector &g) override {
    if (auto *JFNK_mg = dynamic_cast<JFNK_Multigrid<Matrix, Vector> *>(
            this->mf_linear_solver_.get())) {
      return JFNK_mg->initialize(x, g);
    } else {
      QuasiNewtonBase<Vector>::initialize_approximation(x, g);
    }
  }

 private:
  void init_memory(const Layout &layout) {
    s.zeros(layout);
    g.zeros(layout);
    y.zeros(layout);
    g_minus.zeros(layout);

    initialized_ = true;
    layout_ = layout;
  }

  Vector g, s, y, g_minus;
  bool initialized_;
  Layout layout_;
};

}  // namespace utopia
#endif  // UTOPIA_SOLVER_QUASI_NEWTON_HPP
