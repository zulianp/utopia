#ifndef UTOPIA_TR_BOUNDS_GELMAN_MANDEL_HPP
#define UTOPIA_TR_BOUNDS_GELMAN_MANDEL_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_NonLinearSolver.hpp"

#include "utopia_IdentityTransfer.hpp"

#include <iomanip>
#include <limits>

namespace utopia {
template <class Matrix, class Vector>
class TRBoundsGelmanMandel
    : public MultilevelVariableBoundSolverInterface<
          Matrix, Vector, TRBoundsGelmanMandel<Matrix, Vector>> {
 public:
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  typedef utopia::MultilevelVariableBoundSolverInterface<
      Matrix, Vector, TRBoundsGelmanMandel<Matrix, Vector>>
      Base;

  TRBoundsGelmanMandel(
      const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> &transfer)
      : Base(transfer) {}

  void init_memory_impl(const std::vector<Layout> &layouts) {
    constraints_memory_.init_memory(layouts);
  }

  void init_level_impl(const SizeType &level, const Vector &x_finer_level,
                       const Vector &x_level, const Scalar &delta_fine) {
    auto finer_level = level + 1;
    Scalar I_inf_norm = this->transfer_[level]->projection_inf_norm();

    {
      auto d_x_finer = const_local_device_view(x_finer_level);
      auto d_tr_lb = const_local_device_view(
          constraints_memory_.active_lower[finer_level]);
      auto d_help = local_device_view(this->help_[finer_level]);

      parallel_for(local_range_device(this->help_[finer_level]),
                   UTOPIA_LAMBDA(const SizeType i) {
                     auto val = d_x_finer.get(i) - delta_fine;
                     auto lbi = d_tr_lb.get(i);
                     d_help.set(i, device::max(lbi, val));
                   });
    }

    this->help_[finer_level] = this->help_[finer_level] - x_finer_level;
    Scalar lower_multiplier = 1.0 / I_inf_norm * max(this->help_[finer_level]);
    {
      auto d_x = const_local_device_view(x_level);
      auto d_alb = local_device_view(constraints_memory_.active_lower[level]);

      parallel_for(local_range_device(constraints_memory_.active_lower[level]),
                   UTOPIA_LAMBDA(const SizeType i) {
                     d_alb.set(i, d_x.get(i) + lower_multiplier);
                   });
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    {
      auto d_x_finer = const_local_device_view(x_finer_level);
      auto d_tr_ub = const_local_device_view(
          constraints_memory_.active_upper[finer_level]);
      auto d_help = this->help_[finer_level];

      parallel_for(local_range_device(this->help_[finer_level]),
                   UTOPIA_LAMBDA(const SizeType i) {
                     auto val = d_x_finer.get(i) + delta_fine;
                     auto lbi = d_tr_ub.get(i);
                     d_help.set(i, device::min(lbi, val));
                   });
    }

    this->help_[finer_level] = this->help_[finer_level] - x_finer_level;
    Scalar upper_multiplier = 1.0 / I_inf_norm * min(this->help_[finer_level]);

    {
      auto d_x = const_local_device_view(x_level);
      auto d_aub = local_device_view(constraints_memory_.active_upper[level]);

      parallel_for(local_range_device(constraints_memory_.active_upper[level]),
                   UTOPIA_LAMBDA(const SizeType i) {
                     d_aub.set(i, d_x.get(i) + upper_multiplier);
                   });
    }
  }

  const Vector &active_upper(const SizeType &level) {
    return constraints_memory_.active_upper[level];
  }

  const Vector &active_lower(const SizeType &level) {
    return constraints_memory_.active_lower[level];
  }

 private:
  ConstraintsLevelMemory<Vector> constraints_memory_;
};

}  // namespace utopia
#endif  // UTOPIA_TR_BOUNDS_GELMAN_MANDEL_HPP
