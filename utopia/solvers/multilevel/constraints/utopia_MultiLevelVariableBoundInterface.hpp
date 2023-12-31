#ifndef UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
#define UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_LevelMemory.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_MultiLevelVariableBoundInterface.hpp"
#include "utopia_NonLinearSolver.hpp"

#include "utopia_IdentityTransfer.hpp"

#include <iomanip>
#include <limits>

namespace utopia {
template <class Matrix, class Vector, class Child>
class MultilevelVariableBoundSolverInterface {
  using Scalar = typename Traits<Vector>::Scalar;
  using SizeType = typename Traits<Vector>::SizeType;
  using Layout = typename Traits<Vector>::Layout;

  using BoxConstraints = utopia::BoxConstraints<Vector>;

  using Device = typename Traits<Vector>::Device;

public:
  MultilevelVariableBoundSolverInterface(
      const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> &transfer)
      : has_box_constraints_(false), transfer_(transfer) {}

  virtual ~MultilevelVariableBoundSolverInterface() = default;

  virtual void set_box_constraints(BoxConstraints &box) {
    box_constraints_ = box;
    has_box_constraints_ = true;
  }

  virtual BoxConstraints &get_box_constraints() { return box_constraints_; }

  const Vector &active_upper(const SizeType &level) {
    return static_cast<Child *>(this)->active_upper(level);
  }

  const Vector &active_lower(const SizeType &level) {
    return static_cast<Child *>(this)->active_lower(level);
  }

public:
  void init_memory(const std::vector<Layout> &layouts) {
    const SizeType n_levels = layouts.size();
    help_.resize(n_levels);

    for (SizeType l = 0; l < n_levels; l++) {
      help_[l].zeros(layouts[l]);
    }

    return static_cast<Child *>(this)->init_memory_impl(layouts);
  }

  virtual bool check_feasibility(const SizeType &level, Vector &x) {
    SizeType n_terminates = 0;
    // active lower/upper
    {
      auto d_u = const_device_view(active_upper(level));
      auto d_l = const_device_view(active_lower(level));
      auto d_x = const_device_view(x);

      Device::parallel_reduce(
          range(x),
          UTOPIA_LAMBDA(const SizeType i)->SizeType {
            const Scalar xi = d_x.get(i);
            const Scalar li = d_l.get(i);
            const Scalar ui = d_u.get(i);

            return static_cast<SizeType>(xi < li || xi > ui);
          },
          n_terminates);
    }

    bool terminate = n_terminates > 0;
    return x.comm().disjunction(terminate);

    // get_projection(active_lower(level), active_upper(level), x);
    // return false;
  }

  // TODO:: seems redundant with VariableBOundSolverInterface functionality =>
  // remove
  void get_projection(const Vector &lb, const Vector &ub, Vector &x) {
    assert(!empty(lb));
    assert(!empty(ub));
    assert(!empty(x));
    utopia_assert_equal(lb.size(), ub.size());
    utopia_assert_equal(x.size(), ub.size());

    {
      auto d_lb = const_local_view_device(lb);
      auto d_ub = const_local_view_device(ub);
      auto x_view = local_view_device(x);

      parallel_for(
          local_range_device(x), UTOPIA_LAMBDA(const SizeType &i) {
            Scalar li = d_lb.get(i);
            Scalar ui = d_ub.get(i);

            auto xi = x_view.get(i);

            if (li >= xi) {
              x_view.set(i, li);
            } else {
              x_view.set(i, (ui <= xi) ? ui : xi);
            }
          });
    }
  }

  virtual Scalar criticality_measure_inf(const SizeType &level, const Vector &x,
                                         const Vector &g) {
    help_[level] = x - g;
    get_projection(active_lower(level), active_upper(level), help_[level]);
    help_[level] -= x;

    return norm2(help_[level]);
  }

  void init_level(const SizeType &level, const Vector &x_finer_level,
                  const Vector &x_level, const Scalar &delta_fine) {
    return static_cast<Child *>(this)->init_level_impl(level, x_finer_level,
                                                       x_level, delta_fine);
  }

protected:
  BoxConstraints box_constraints_; // constraints on the finest level....
  bool has_box_constraints_; // as we can run rmtr with inf. norm also without
                             // constraints...

  const std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> &transfer_;

protected:
  std::vector<Vector> help_;
};

} // namespace utopia
#endif // UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
