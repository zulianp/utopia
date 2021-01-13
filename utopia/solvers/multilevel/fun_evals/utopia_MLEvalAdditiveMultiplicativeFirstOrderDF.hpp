#ifndef UTOPIA_ML_EVAL_ADDITIVE_FIRST_ORDER_MULTIPLICATIVE_DERIVATIVE_FREE_HPP
#define UTOPIA_ML_EVAL_ADDITIVE_FIRST_ORDER_MULTIPLICATIVE_DERIVATIVE_FREE_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia {
// Matrix free first order
template <typename Matrix, typename Vector>
class MultilevelDerivEval<Matrix, Vector,
                          FIRST_ORDER_ADDITIVE_MULTIPLICATIVE_DF>
    final {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

 public:
  MultilevelDerivEval(const SizeType &nl_levels)
      : n_levels_(nl_levels), initialized_(false), gamma_(0.5) {}

  void gamma(const Scalar &gamma) { gamma_ = gamma; }

  Scalar gamma() { return gamma_; }

  inline Scalar compute_energy(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global) {
    Scalar energy = 0.0;
    Scalar energy_m = 0.0;
    Scalar energy_a = 0.0;

    fun.value(x, energy);

    if (level < n_levels_ - 1) {
      energy_m = (e_diff_m[level] + dot(g_diff_m[level], s_global)) * energy;
      energy_a = energy + e_diff_a[level] + dot(g_diff_a[level], s_global);

      energy = gamma_ * energy_a + ((1.0 - gamma_) * energy_m);
    }

    return energy;
  }

  // s_global is assummed to be zero
  inline Scalar compute_energy(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x) {
    Scalar energy = 0.0;
    Scalar energy_m = 0.0;
    Scalar energy_a = 0.0;
    fun.value(x, energy);
    if (level < n_levels_ - 1) {
      energy_m = (e_diff_m[level]) * energy;
      energy_a = energy + e_diff_a[level];

      energy = gamma_ * energy_a + ((1.0 - gamma_) * energy_m);
    }
    return energy;
  }

  inline bool compute_gradient(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global) {
    fun.gradient(x, g[level]);
    Scalar energy = 0.0;
    fun.value(x, energy);

    if (level < n_levels_ - 1) {
      gm[level] = g[level] * (e_diff_m[level] + dot(g_diff_m[level], s_global));
      gm[level] += g_diff_m[level] * energy;

      ga[level] = g[level] + g_diff_a[level];
      g[level] = gamma_ * ga[level] + ((1.0 - gamma_) * gm[level]);
    }

    return true;
  }

  // s_global is assummed to be zero
  inline bool compute_gradient(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x) {
    fun.gradient(x, g[level]);
    Scalar energy = 0.0;
    fun.value(x, energy);

    if (level < n_levels_ - 1) {
      gm[level] = g[level] * e_diff_m[level];
      gm[level] += g_diff_m[level] * energy;

      ga[level] = g[level] + g_diff_a[level];
      g[level] = gamma_ * ga[level] + ((1.0 - gamma_) * gm[level]);
    }

    return true;
  }

  inline bool compute_gradient(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global,
                               const Scalar &energy) {
    fun.gradient(x, g[level]);

    if (level < n_levels_ - 1) {
      gm[level] = g[level] * e_diff_m[level];
      gm[level] += g_diff_m[level] * energy;

      ga[level] = g[level] + g_diff_a[level];
      g[level] = gamma_ * ga[level] + ((1.0 - gamma_) * gm[level]);
    }

    return true;
  }

  inline Scalar compute_gradient_energy(
      const SizeType &level, const ExtendedFunction<Matrix, Vector> &fun,
      const Vector &x, const Vector &s_global) {
    Scalar energy = 0.0;
    Scalar energy_m = 0.0;
    Scalar energy_a = 0.0;
    fun.value(x, energy);
    fun.gradient(x, g[level]);

    if (level < n_levels_ - 1) {
      energy_m = (e_diff_m[level] + dot(g_diff_m[level], s_global)) * energy;
      energy_a = energy + e_diff_a[level] + dot(g_diff_a[level], s_global);

      energy = gamma_ * energy_a + ((1.0 - gamma_) * energy_m);

      gm[level] = g[level] * (e_diff_m[level] + dot(g_diff_m[level], s_global));
      gm[level] += g_diff_m[level] * energy;

      ga[level] = g[level] + g_diff_a[level];
      g[level] = gamma_ * ga[level] + ((1.0 - gamma_) * gm[level]);
    }

    return energy;
  }

  void init_memory(
      const std::vector<Layout> &layouts,
      const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > >
          & /*level_functions*/) {
    g_diff_m.resize(n_levels_);
    g_diff_a.resize(n_levels_);
    g.resize(n_levels_);
    ga.resize(n_levels_);
    gm.resize(n_levels_);
    y.resize(n_levels_);
    e_diff_m.resize(n_levels_);
    e_diff_a.resize(n_levels_);

    for (auto l = 0; l < n_levels_; l++) {
      g_diff_m[l].zeros(layouts[l]);
      g_diff_a[l].zeros(layouts[l]);
      g[l].zeros(layouts[l]);
      ga[l].zeros(layouts[l]);
      gm[l].zeros(layouts[l]);
      y[l].zeros(layouts[l]);
      e_diff_m[l] = 0.0;
      e_diff_a[l] = 0.0;
    }

    initialized_ = true;
  }

  bool initialized() const { return initialized_; }

 private:
  SizeType n_levels_;
  bool initialized_;
  Scalar gamma_;

 public:
  std::vector<Vector> g, ga, gm, g_diff_m, y, g_diff_a;
  std::vector<Scalar> e_diff_m, e_diff_a;
};
}  // namespace utopia

#endif  // UTOPIA_ML_EVAL_FIRST_ORDER_DERIV_FREE_HPP
