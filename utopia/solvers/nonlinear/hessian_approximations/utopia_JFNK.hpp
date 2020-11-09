#ifndef UTOPIA_JFNK_HPP
#define UTOPIA_JFNK_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"

namespace utopia {

template <class Vector>
class JFNK final : public HessianApproximation<Vector> {
  using Scalar = typename utopia::Traits<Vector>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

 public:
  JFNK(const FunctionBase<Vector> &fun) : eps_(1e-14), fun_(fun) {}

  void initialize(const Vector &x_k, const Vector &g) override {
    HessianApproximation<Vector>::initialize(x_k, g);

    x_k_ = x_k;
    g_ = g;

    ones_ = 0.0 * g_;
    ones_.set(1.0);

    this->initialized(true);
  }

  void reset() override {
    Vector x, g;
    this->initialize(x, g);
  }

  inline JFNK<Vector> *clone() const override {
    return new JFNK<Vector>(*this);
  }

  bool update(const Vector & /*s*/, const Vector & /*y*/, const Vector &x,
              const Vector &g) override {
    x_k_ = x;
    g_ = g;

    return true;
  }

  bool apply_Hinv(const Vector & /*g*/, Vector & /*q*/) override {
    utopia_error("utopia::JFNK::apply_Hinv:: not supported... \n");
    return false;
  }

  bool apply_H(const Vector &v, Vector &result) override {
    Vector aa = std::sqrt(eps_) * (ones_ + x_k_);
    Scalar sum_a = sum(aa);

    SizeType n_glob = size(result).get(0);

    Scalar per_ = 0.0;

    if (norm2(v) > eps_) {
      per_ = 1. / (n_glob * Scalar(norm2(v)));
      per_ *= sum_a;
    } else {
      per_ = sum_a / n_glob;
    }

    Vector x_p = x_k_ + (per_ * v);

    Vector grad_pertubed;
    fun_.gradient(x_p, grad_pertubed);

    result = (1. / per_) * (grad_pertubed - g_);

    return true;
  }

  void read(Input &in) override { HessianApproximation<Vector>::read(in); }

  void print_usage(std::ostream &os) const override {
    HessianApproximation<Vector>::print_usage(os);
  }

 private:
  Scalar eps_;
  const FunctionBase<Vector> &fun_;
  Vector x_k_, g_, ones_;
};

}  // namespace utopia

#endif  // UTOPIA_JFNK_HPP
