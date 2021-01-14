#ifndef UTOPIA_EXTENDED_FUNCTION_HPP
#define UTOPIA_EXTENDED_FUNCTION_HPP

#include <utility>

#include "utopia_Base.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
/**
 * @brief      Class for Nonlinear Function, all application context needed by
 * solver is usually provided inside of this functions. In optimization
 * settings, user needs to supply value(energy), gradient, hessian. Difference,
 * between Function and ExtendedFunction is that here, we make use of additional
 * informations to improve convergence
 *
 * @tparam     Matrix
 * @tparam     Vector
 */
template <class Matrix, class Vector>
class ExtendedFunction : public Function<Matrix, Vector> {
 public:
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  ~ExtendedFunction() override = default;

  ExtendedFunction() = default;

  ExtendedFunction(const Vector &x_init, const Vector &bc_marker) {
    this->set_equality_constrains(bc_marker, x_init);
  }

  bool value(const Vector & /*point*/, Scalar & /*value*/) const override = 0;

  // Copy of vec...
  virtual Vector initial_guess() const { return _x_eq_values; }

  virtual Layout layout() const { return utopia::layout(_x_eq_values); }

  bool hessian(const Vector &x, Matrix &H) const override = 0;
  bool hessian(const Vector & /*point*/, Matrix & /*result*/,
               Matrix & /*preconditioner*/) const override {
    return false;
  }

  bool has_preconditioner() const override { return false; }

  bool update(const Vector & /*point*/) override { return true; }

  virtual bool get_eq_constrains_values(Vector &x) const {
    x = _x_eq_values;
    return true;
  }

  virtual bool get_eq_constrains_flg(Vector &x) const {
    x = _eq_constrains_flg;
    return true;
  }

  inline const std::vector<SizeType> &get_indices_related_to_BC() const {
    return indices_eq_constraints_;
  }

  Vector &get_eq_constrains_flg() { return _eq_constrains_flg; }

  Vector &get_eq_constrains_values() { return _x_eq_values; }

  virtual bool set_equality_constrains(const Vector &eq_constrains_flg,
                                       const Vector &x_in) {
    _x_eq_values = x_in;
    _eq_constrains_flg = eq_constrains_flg;

    this->init_constraint_indices();

    return true;
  }

  bool init_constraint_indices() {
    indices_eq_constraints_.clear();

    {
      Read<Vector> r(_eq_constrains_flg);

      Range range_w = range(_eq_constrains_flg);
      for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
        if (_eq_constrains_flg.get(i) == 1) {
          indices_eq_constraints_.push_back(i);
        }
      }
    }

    return true;
  }

  virtual bool zero_contribution_to_equality_constrains(Vector &x) const {
    UTOPIA_NO_ALLOC_BEGIN("RMTR::zero_contribution_to_equality_constrains");
    // x = _eq_constraints_mask_matrix_ * x;

    {
      auto d_flg = const_local_view_device(_eq_constrains_flg);
      auto x_view = local_view_device(x);

      parallel_for(local_range_device(x), UTOPIA_LAMBDA(const SizeType &i) {
        Scalar flg = d_flg.get(i);

        // TODO:: use abs with eps tolerance
        if (flg == 1.0) {
          x_view.set(i, 0.0);
        }
      });
    }

    UTOPIA_NO_ALLOC_END();
    return true;
  }

 protected:
  Vector _x_eq_values;
  Vector _eq_constrains_flg;

  std::vector<SizeType> indices_eq_constraints_;
};

template <class Matrix, class Vector>
class Function_rhs final : public ExtendedFunction<Matrix, Vector> {
 public:
  using Scalar = typename utopia::Traits<Matrix>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

  typedef utopia::ExtendedFunction<Matrix, Vector> Fun;

  Function_rhs(std::shared_ptr<Fun> fun) : fun_(std::move(fun)) {}

  bool value(const Vector &x, Scalar &value) const override {
    fun_->value(x, value);
    return true;
  }

  bool gradient(const Vector &x, Vector &g) const override {
    fun_->gradient(x, g);

    if (size(g) == size(this->rhs_)) {
      g = g - this->rhs_;
    }

    return true;
  }

  bool hessian(const Vector &x, Matrix &H) const override {
    fun_->hessian(x, H);
    return true;
  }

  virtual bool set_rhs(const Vector &rhs) {
    rhs_ = rhs;
    return true;
  }

  virtual bool reset_rhs() {
    if (!empty(rhs_)) {
      rhs_.set(0.0);
    } else {
      utopia_error("error in reset rhs... \n");
    }
    return true;
  }

  virtual bool get_rhs(Vector &rhs) const {
    rhs = rhs_;
    return true;
  }

  virtual bool has_rhs() const { return !empty(rhs_); }

 private:
  std::shared_ptr<Fun> fun_;
  Vector rhs_;
};

}  // namespace utopia
#endif  // UTOPIA_EXTENDED_FUNCTION_HPP
