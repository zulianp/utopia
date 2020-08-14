#ifndef UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP
#define UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia {
template <typename Matrix, typename Vector>
class MultilevelDerivEval<Matrix, Vector, SECOND_ORDER_DF> final {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  using HessianApprox = utopia::HessianApproximation<Vector>;
  using Operator = utopia::Operator<Vector>;

 public:
  MultilevelDerivEval(const SizeType &nl_levels)
      : n_levels_(nl_levels), initialized_(false) {
    utopia::out() << "this class works only with the identity transfer \n";
  }

  inline Scalar compute_energy(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global) {
    Scalar energy = 0.0;
    fun.value(x, energy);

    if (level < n_levels_ - 1) {
      // help_[level] = H_diff[level] * s_global;
      hessian_approxs_init_[level][0]->apply_H(s_global, help1_[level]);
      hessian_approxs_init_[level][1]->apply_H(s_global, help2_[level]);

      help1_[level] -= help2_[level];

      energy +=
          (0.5 * dot(help1_[level], s_global)) + dot(g_diff[level], s_global);
    }
    return energy;
  }

  // s_global is assummed to be zero
  inline Scalar compute_energy(const SizeType & /*level*/,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x) {
    Scalar energy = 0.0;
    fun.value(x, energy);
    return energy;
  }

  inline bool compute_gradient(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global) {
    fun.gradient(x, g[level]);

    if (level < n_levels_ - 1) {
      // help_[level] = H_diff[level] * s_global;

      hessian_approxs_init_[level][0]->apply_H(s_global, help1_[level]);
      hessian_approxs_init_[level][1]->apply_H(s_global, help2_[level]);

      help1_[level] -= help2_[level];

      g[level] = g[level] + g_diff[level] + help1_[level];
    }
    return true;
  }

  // s_global is assummed to be zero
  inline bool compute_gradient(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x) {
    fun.gradient(x, g[level]);

    if (level < n_levels_ - 1) {
      g[level] += g_diff[level];
    }
    return true;
  }

  inline Scalar compute_gradient_energy(
      const SizeType &level, const ExtendedFunction<Matrix, Vector> &fun,
      const Vector &x, const Vector &s_global) {
    Scalar energy = 0.0;
    fun.value(x, energy);
    fun.gradient(x, g[level]);

    if (level < n_levels_ - 1) {
      // help_[level] = H_diff[level] * s_global;

      hessian_approxs_init_[level][0]->apply_H(s_global, help1_[level]);
      hessian_approxs_init_[level][1]->apply_H(s_global, help2_[level]);

      help1_[level] -= help2_[level];

      energy +=
          (0.5 * dot(help1_[level], s_global)) + dot(g_diff[level], s_global);
      g[level] += g_diff[level] + help1_[level];
    }

    return energy;
  }

  inline bool compute_hessian(const SizeType &level,
                              const ExtendedFunction<Matrix, Vector> &fun,
                              const Vector &x) {
    // fun.hessian(x, H[level]);

    // if (level < n_levels_ - 1) {
    //   H[level] += H_diff[level];
    // }

    utopia::err()
        << "compute hessian should not be used for matrix free approach ";

    return true;
  }

  std::shared_ptr<Operator> build_apply_H_diff() {
    std::function<void(const SizeType &, const Vector &, Vector &)> my_func =
        [this](const SizeType &level, const Vector &x, Vector &result) {
          // this->apply_H(x, result);

          hessian_approxs_init_[level][0]->apply_H(x, help1_[level]);
          hessian_approxs_init_[level][1]->apply_H(x, help2_[level]);

          result = help1_[level] - help2_[level];
        };

    return std::make_shared<Operator>(*this, my_func);
  }

  std::shared_ptr<Operator> build_apply_H_plus_Hdiff(
      const std::shared_ptr<Operator> &apply_B_fun) {
    std::function<void(const SizeType &, const Vector &, Vector &)> my_func =
        [this, &apply_B_fun](const SizeType &level, const Vector &x,
                             Vector &result) {
          // this->apply_H(x, result);

          hessian_approxs_init_[level][0]->apply_H(x, help1_[level]);
          hessian_approxs_init_[level][1]->apply_H(x, help2_[level]);

          result = help1_[level] - help2_[level];

          apply_B_fun->apply(help1_[level]);
          result += help1_[level];

        };

    return std::make_shared<Operator>(*this, my_func);
  }

  void init_memory(
      const std::vector<Layout> &layouts,
      const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector>>>
          &level_functions) {
    help1_.resize(n_levels_);
    help2_.resize(n_levels_);
    g_diff.resize(n_levels_);
    g.resize(n_levels_);

    // this one could be done nicer
    hessian_approxs_init_.resize(n_levels_);

    for (auto l = 0; l < n_levels_; l++) {
      help1_[l].zeros(layouts[l]);
      help2_[l].zeros(layouts[l]);
      g_diff[l].zeros(layouts[l]);
      g[l].zeros(layouts[l]);

      hessian_approxs_init_[l].resize(2);
    }

    initialized_ = true;
  }

  bool initialized() const { return initialized_; }

  void set_init_approx(const SizeType &level, const HessianApprox &approx_fine,
                       const HessianApprox &approx_coarse) {
    hessian_approxs_init_[level, 0] = approx_fine;
    hessian_approxs_init_[level, 1] = approx_coarse;
  }

 private:
  SizeType n_levels_;
  bool initialized_;

 public:
  std::vector<Vector> g, g_diff, help1_, help2_;

  //          fine_presmoothing, coarse_init
  std::vector<std::vector<HessianApprox>> hessian_approxs_init_;
};

}  // namespace utopia

#endif  // UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP