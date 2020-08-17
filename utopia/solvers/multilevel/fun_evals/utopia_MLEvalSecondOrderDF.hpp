#ifndef UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP
#define UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_HessianApproximations.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

// TODO:: extend for generic transfer
// TODO:: use help1_, help2_ based on level informations
namespace utopia {
template <typename Matrix, typename Vector>
class MultilevelDerivEval<Matrix, Vector, SECOND_ORDER_DF> final {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;
  using Communicator = typename Traits<Vector>::Communicator;

  using HessianApprox = utopia::HessianApproximation<Vector>;
  using HessianApproxPtr = std::shared_ptr<HessianApprox>;
  using FunctionOperatorHA = typename HessianApprox::FunctionOperator;

  using Operator = utopia::Operator<Vector>;

 public:
  struct HessianApproxPtrConsistency {
   public:
    HessianApproxPtrConsistency(const HessianApproxPtr &fine,
                                const HessianApproxPtr &coarse) {
      fine_ = std::shared_ptr<HessianApprox>(fine->clone());
      coarse_ = std::shared_ptr<HessianApprox>(coarse->clone());
    }

    HessianApproxPtrConsistency() {}

    HessianApproxPtr fine_;
    HessianApproxPtr coarse_;
  };

  class FunctionOperatorDFEval final : public Operator {
   public:
    FunctionOperatorDFEval(
        MultilevelDerivEval &parent, const Size size, const Size local_size,
        const std::shared_ptr<Communicator> comm,
        const std::function<void(const Vector &, Vector &)> &operator_action)
        : parent_(parent),
          size_(size),
          local_size_(local_size),
          comm_(comm),
          operator_action_(operator_action) {}

    bool apply(const Vector &rhs, Vector &ret) const override {
      operator_action_(rhs, ret);
      return true;
    }

    inline Communicator &comm() override { return *comm_; }

    inline const Communicator &comm() const override { return *comm_; }

    inline Size size() const override { return size_; }

    inline Size local_size() const override { return local_size_; }

   private:
    MultilevelDerivEval &parent_;
    Size size_, local_size_;
    std::shared_ptr<Communicator> comm_;
    const std::function<void(const Vector &, Vector &)> operator_action_;
  };

  MultilevelDerivEval(const SizeType &nl_levels)
      : n_levels_(nl_levels), initialized_(false) {
    utopia::out() << "this class works only with the identity transfer \n";
  }

  ~MultilevelDerivEval() = default;

  inline Scalar compute_energy(const SizeType &level,
                               const ExtendedFunction<Matrix, Vector> &fun,
                               const Vector &x, const Vector &s_global) {
    Scalar energy = 0.0;
    fun.value(x, energy);

    if (level < n_levels_ - 1) {
      // help_[level] = H_diff[level] * s_global;
      hessian_approxs_init_[level].fine_->apply_H(x, help1_[level]);
      hessian_approxs_init_[level].coarse_->apply_H(x, help2_[level]);

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

      hessian_approxs_init_[level].fine_->apply_H(x, help1_[level]);
      hessian_approxs_init_[level].coarse_->apply_H(x, help2_[level]);

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

      hessian_approxs_init_[level].fine_->apply_H(x, help1_[level]);
      hessian_approxs_init_[level].coarse_->apply_H(x, help2_[level]);

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

  std::shared_ptr<FunctionOperatorDFEval> build_apply_H_diff(
      const SizeType &level) {
    std::function<void(const Vector &, Vector &)> my_func =
        [this, level](const Vector &x, Vector &result) {

          // std::cout << "------ inside1 ------- \n";
          this->hessian_approxs_init_[level].fine_->apply_H(x, help1_[level]);

          // std::cout << "------ inside2 ------- \n";
          this->hessian_approxs_init_[level].coarse_->apply_H(x, help2_[level]);

          result = this->help1_[level] - this->help2_[level];
          // std::cout << "------ inside2 ------- \n";

        };

    auto comm = std::shared_ptr<Communicator>(help1_[level].comm().clone());

    Size size, local_size;
    size.set_dims(2);
    size.set(0, help1_[level].size());
    size.set(1, help1_[level].size());

    local_size.set_dims(2);
    local_size.set(0, help1_[level].local_size());
    local_size.set(1, help1_[level].local_size());

    return std::make_shared<FunctionOperatorDFEval>(*this, size, local_size,
                                                    comm, my_func);
  }

  std::shared_ptr<FunctionOperatorDFEval> build_apply_H_plus_Hdiff(
      const SizeType &level,
      const std::shared_ptr<FunctionOperatorHA> &apply_B_fun) {
    std::function<void(const Vector &, Vector &)> my_func =
        [this, &apply_B_fun, level](const Vector &x, Vector &result) {
          // this->apply_H(x, result);

          // std::cout << "fine:  \n";
          this->hessian_approxs_init_[level].fine_->apply_H(x, help1_[level]);

          // std::cout << "coarse:  \n";
          this->hessian_approxs_init_[level].coarse_->apply_H(x, help2_[level]);

          // std::cout << "sum:  \n";
          result = this->help1_[level] - this->help2_[level];

          // std::cout << "new-one:  \n";
          apply_B_fun->apply(x, help1_[level]);
          result += this->help1_[level];

          // std::cout << "--------------------------- \n";
        };

    auto comm = std::shared_ptr<Communicator>(help1_[level].comm().clone());

    Size size, local_size;
    size.set_dims(2);
    size.set(0, help1_[level].size());
    size.set(1, help1_[level].size());

    local_size.set_dims(2);
    local_size.set(0, help1_[level].local_size());
    local_size.set(1, help1_[level].local_size());

    return std::make_shared<FunctionOperatorDFEval>(*this, size, local_size,
                                                    comm, my_func);
  }

  void init_memory(
      const std::vector<Layout> &layouts,
      const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector>>>
          &level_functions) {
    help1_.resize(n_levels_);
    help2_.resize(n_levels_);
    g_diff.resize(n_levels_);
    g.resize(n_levels_);
    y.resize(n_levels_);

    // this one could be done nicer
    hessian_approxs_init_.resize(n_levels_);

    for (auto l = 0; l < n_levels_; l++) {
      help1_[l].zeros(layouts[l]);
      help2_[l].zeros(layouts[l]);
      g_diff[l].zeros(layouts[l]);
      g[l].zeros(layouts[l]);
      y[l].zeros(layouts[l]);
    }

    initialized_ = true;
  }

  bool initialized() const { return initialized_; }

  void set_init_approx(const SizeType &level,
                       const HessianApproxPtr &approx_fine,
                       const HessianApproxPtr &approx_coarse) {
    // does not need to be cloned
    hessian_approxs_init_[level].fine_ =
        std::shared_ptr<HessianApprox>(approx_fine->clone());

    // we need to clone as we want to stop updates
    hessian_approxs_init_[level].coarse_ =
        std::shared_ptr<HessianApprox>(approx_coarse->clone());
  }

 private:
  SizeType n_levels_;
  bool initialized_;

 public:
  std::vector<Vector> g, g_diff, help1_, help2_, y;

  //          fine_presmoothing, coarse_init
  std::vector<HessianApproxPtrConsistency> hessian_approxs_init_;
};

}  // namespace utopia

#endif  // UTOPIA_ML_EVAL_SECOND_ORDER_DF_HPP