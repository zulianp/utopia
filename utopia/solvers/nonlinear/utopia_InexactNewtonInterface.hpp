#ifndef UTOPIA_NEWTON_INEXACT_FORCING_NEWTON_INTERFACE_HPP
#define UTOPIA_NEWTON_INEXACT_FORCING_NEWTON_INTERFACE_HPP

#include <utility>

#include "utopia_ConjugateGradient.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_Function.hpp"
#include "utopia_Monitor.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_PrintInfo.hpp"

namespace utopia {

enum InexactNewtonForcingStartegies {
  ZERO = 0,
  CAI = 1,
  DEMBO = 2,
  SUPERLINEAR = 3,
  QUADRATIC = 4,
  QUADRATIC_2 = 5
};

template <class Vector>
class InexactNewtonInterface : virtual public Configurable {
 public:
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;

  InexactNewtonInterface()
      : forcing_strategy_(InexactNewtonForcingStartegies::ZERO) {}

  ~InexactNewtonInterface() override = default;

  virtual void forcing_strategy(
      const InexactNewtonForcingStartegies &strategy) {
    forcing_strategy_ = strategy;
  }

  virtual InexactNewtonForcingStartegies forcing_strategy() {
    return forcing_strategy_;
  }

  virtual bool has_forcing_strategy() {
    return (forcing_strategy_ > 0) ? true : false;
  }

 protected:
  virtual Scalar estimate_ls_atol(const Scalar &gnorm, const Scalar &it) {
    Scalar atol_new = 1e-11;

    if (forcing_strategy_ == InexactNewtonForcingStartegies::CAI) {
      atol_new = 10e-4 * gnorm;
    } else if (forcing_strategy_ == InexactNewtonForcingStartegies::DEMBO) {
      atol_new = (gnorm * std::min(std::sqrt(gnorm), 1. / (it + 1)));
    } else if (forcing_strategy_ ==
               InexactNewtonForcingStartegies::SUPERLINEAR) {
      atol_new = (gnorm * std::min(std::sqrt(gnorm), 0.5));
    } else if (forcing_strategy_ == InexactNewtonForcingStartegies::QUADRATIC) {
      atol_new = (gnorm * std::min(gnorm, 1. / (it + 1)));
    } else if (forcing_strategy_ ==
               InexactNewtonForcingStartegies::QUADRATIC_2) {
      atol_new = (gnorm * std::min(gnorm, 0.5));
    } else {
      utopia_error(
          "utopia::Newton:: invalid choice of forcing strategy.... \n");
    }

    atol_new = atol_new < 1e-11 ? 1e-11 : atol_new;

    return atol_new;
  }

  InexactNewtonForcingStartegies forcing_strategy_;
};
}  // namespace utopia

#endif  // UTOPIA_NEWTON_INEXACT_FORCING_NEWTON_INTERFACE_HPP
