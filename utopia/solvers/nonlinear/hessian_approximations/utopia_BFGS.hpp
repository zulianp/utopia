#ifndef UTOPIA_HESSIAN_BFGS_HPP
#define UTOPIA_HESSIAN_BFGS_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {

template <class Matrix, class Vector>
class BFGS final : public HessianApproximation<Vector> {
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

 public:
  BFGS() : HessianApproximation<Vector>(), current_it_(0) {}

  inline BFGS<Matrix, Vector> *clone() const override {
    return new BFGS<Matrix, Vector>(*this);
  }

  void initialize(const Vector & /* x */, const Vector & /* g */) override {
    current_it_ = 0;
    this->initialized(true);

    H_prev_ = Matrix();
    H_prev_inv_ = Matrix();
  }

  void reset() override {
    if (!empty(H_prev_)) H_prev_.identity(layout(H_prev_));

    if (!empty(H_prev_inv_)) H_prev_inv_.identity(layout(H_prev_inv_));
  }

  bool update(const Vector &s_in, const Vector &y_in, const Vector & /* x */,
              const Vector & /* g */) override {
    if (!this->initialized()) {
      utopia_error(
          "BFGS::update: Initialization needs to be done before updating. \n");
      return false;
    }

    if (current_it_ == 0) {
      // SizeType n = local_size(s_in).get(0);

      auto mat_layout = square_matrix_layout(layout(s_in));

      if (update_hessian_) {
        H_prev_.identity(mat_layout, 1.0);
      }

      H_prev_inv_.identity(mat_layout, 1.0);
    }

    s_ = s_in;
    y_ = y_in;

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:1");
    this->update_Hessian_inverse();
    UTOPIA_NO_ALLOC_END();

    if (update_hessian_) {
      UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:2");
      this->update_Hessian();
      UTOPIA_NO_ALLOC_END();
    }

    current_it_++;

    return true;
  }

  Matrix &get_Hessian() { return H_prev_; }

  bool apply_Hinv(const Vector &g, Vector &s) override {
    if (this->initialized()) {
      if (!empty(H_prev_inv_))
        s = H_prev_inv_ * g;
      else
        s = g;
    } else {
      utopia_error(
          "BFGS::apply_Hinv: Initialization needs to be done first. \n");
    }

    return true;
  }

  bool apply_H(const Vector &v, Vector &result) override {
    if (update_hessian_ && this->initialized()) {
      if (!empty(H_prev_))
        result = H_prev_ * v;
      else
        result = v;
    } else {
      utopia_error(
          "BFGS::apply_H can be used only, if H is computed. \n Please turn on "
          "update_hessian option. \n");
    }

    return true;
  }

  Scalar compute_uHinvv_dot(const Vector &u, const Vector &v) override {
    if (this->initialized()) {
      if (!empty(H_prev_inv_))
        return dot(u, H_prev_inv_ * v);
      else
        return dot(u, v);
    } else {
      utopia_error(
          "BFGS::compute_uHinvv_dot: Initialization needs to be done first. "
          "\n");
      return false;
    }
  }

  Scalar compute_uHv_dot(const Vector &u, const Vector &v) override {
    if (update_hessian_ && this->initialized()) {
      if (!empty(H_prev_))
        return dot(u, H_prev_ * v);
      else
        return dot(u, v);
    } else {
      utopia_error(
          "BFGS::compute_uHv_dot can be used only, if H is computed. \n Please "
          "turn on update_hessian "
          "option. \n");
      return 0;
    }
  }

  Scalar compute_uHu_dot(const Vector &u) override {
    if (update_hessian_ && this->initialized()) {
      if (!empty(H_prev_))
        return dot(u, H_prev_ * u);
      else
        return dot(u, u);
    } else {
      utopia_error(
          "BFGS::compute_uHu_dot can be used only, if H is computed. \n Please "
          "turn on update_hessian "
          "option. \n");
      return 0;
    }
  }

  void update_hessian(const bool flg) { update_hessian_ = flg; }

  bool update_hessian() const { return update_hessian_; }

  void read(Input &in) override {
    HessianApproximation<Vector>::read(in);
    in.get("update_hessian", update_hessian_);
  }

  void print_usage(std::ostream &os) const override {
    HessianApproximation<Vector>::print_usage(os);
    this->print_param_usage(os, "update_hessian", "bool", "Default step size.",
                            "false");
  }

 private:
  void update_Hessian_inverse() {
    Scalar s_y = dot(s_, y_);
    Scalar yHy = dot(y_, H_prev_inv_ * y_);

    // checking for numerical instabilities
    if (s_y < this->num_tol() || yHy < this->num_tol() || !std::isfinite(s_y) ||
        !std::isfinite(yHy)) {
      // std::cout<<"--- BFGS::update_Hessian_inverse is reaching num. tolerance
      // \n";
      return;
    }

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:3");
    ss_ = outer(s_, s_);
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:4");
    Hy_ = H_prev_inv_ * y_;
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:5");
    Hys_ = outer(Hy_, s_);
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:6");
    sy_outerH_ = outer(s_, y_);
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:7");
    sy_outerH_ = sy_outerH_ * H_prev_inv_;
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:8");
    H_prev_inv_ += (s_y + yHy) / (s_y * s_y) * ss_;
    UTOPIA_NO_ALLOC_END();

    UTOPIA_NO_ALLOC_BEGIN("Quasi BFGS:9");
    H_prev_inv_ -= (1. / s_y) * (Hys_ + sy_outerH_);
    UTOPIA_NO_ALLOC_END();
  }

  void update_Hessian() {
    Scalar y_s = dot(y_, s_);

    Vector H_s = H_prev_ * s_;
    Scalar sHs = dot(s_, H_s);

    if (y_s < this->num_tol() || !std::isfinite(y_s) || sHs < this->num_tol() ||
        !std::isfinite(sHs)) {
      // std::cout<<"--- BFGS::update_Hessian is reaching num. tolerance \n";
      return;
    }

    Matrix yy = outer(y_, y_);
    H_prev_ += 1. / y_s * yy;

    Matrix a = outer(H_s, H_s);
    H_prev_ -= 1. / sHs * a;
  }

 private:
  Vector s_;  // x_{k+1} - x_{k}
  Vector y_;  // g_{k+1} - g_{k}

  Matrix H_prev_inv_;  // H^{-1}_{k}
  Matrix H_prev_;      // H^{1}_{k}

  bool update_hessian_{false};
  SizeType current_it_;

  // help mats/vecs
  Matrix ss_, Hys_, sy_outerH_;
  Vector Hy_;
};

}  // namespace utopia

#endif  // UTOPIA_HESSIAN_BFGS_HPP