#ifndef UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
#define UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS

#include <iomanip>
#include <locale>
#include "utopia_Base.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"

namespace utopia {

// This function is used for implementation of test functions for unconstrained
// nonlinear benchmark
template <class Matrix, class Vector>
class TestFunctionInterface {
 public:
  using Traits = utopia::Traits<Vector>;
  using Scalar = typename Traits::Scalar;
  using SizeType = typename Traits::SizeType;
  using Comm = typename Traits::Communicator;

  virtual ~TestFunctionInterface() = default;

  virtual Vector initial_guess() const = 0;
  virtual const Vector &exact_sol() const = 0;
  virtual Scalar min_function_value() const = 0;

  virtual std::string name() const = 0;
  virtual SizeType dim() const = 0;

  virtual bool exact_sol_known() const { return true; }

  virtual bool parallel() const { return false; }

  virtual void describe() const {
    if (mpi_world_rank() == 0) {
      utopia::out() << name() << ",   Globalsize:";

      std::string numWithCommas = std::to_string(dim());
      int insertPosition = numWithCommas.length() - 3;
      while (insertPosition > 0) {
        numWithCommas.insert(insertPosition, ",");
        insertPosition -= 3;
      }

      utopia::out() << numWithCommas << ",   parallel:  " << parallel()
                    << ",   sol. known:" << exact_sol_known() << "  \n";
    }
  }
};

template <class Matrix, class Vector>
class UnconstrainedTestFunction
    : virtual public TestFunctionInterface<Matrix, Vector>,
      virtual public Function<Matrix, Vector> {
 public:
  using Scalar = typename utopia::Traits<Matrix>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

  ~UnconstrainedTestFunction() override = default;
};

template <class Matrix, class Vector>
class UnconstrainedExtendedTestFunction
    : virtual public TestFunctionInterface<Matrix, Vector>,
      virtual public ExtendedFunction<Matrix, Vector> {
 public:
  using Scalar = typename utopia::Traits<Matrix>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;

  ~UnconstrainedExtendedTestFunction() override = default;

  Vector initial_guess() const override = 0;
};

template <class Matrix, class Vector>
class ConstrainedTestFunction
    : virtual public TestFunctionInterface<Matrix, Vector>,
      virtual public Function<Matrix, Vector> {
 public:
  using Scalar = typename utopia::Traits<Matrix>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  using BoxConstraints = utopia::BoxConstraints<Vector>;

  ~ConstrainedTestFunction() override = default;

  virtual const BoxConstraints &box_constraints() const { return constraints_; }

  virtual void set_box_constraints(const BoxConstraints &box) {
    constraints_ = box;
  }

  virtual const Vector &upper_bound() const {
    if (!constraints_.upper_bound()) {
      utopia_error("ConstrainedTestFunction::upper bound does not exist. \n");
    }

    return *constraints_.upper_bound();
  }

  virtual const Vector &lower_bound() const {
    if (!constraints_.lower_bound()) {
      utopia_error("ConstrainedTestFunction::lower bound does not exist. \n");
    }

    return *constraints_.lower_bound();
  }

  virtual bool has_lower_bound() const {
    return constraints_.has_lower_bound();
  }

  virtual bool has_upper_bound() const {
    return constraints_.has_upper_bound();
  }

  virtual bool is_feasible(Vector &x) {
    if (!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
      return true;

    if (empty(help_) || size(help_) != size(x)) {
      help_ = 0.0 * x;
    }

    if (constraints_.has_upper_bound() && constraints_.has_lower_bound()) {
      const auto &ub = *constraints_.upper_bound();
      const auto &lb = *constraints_.lower_bound();

      {
        auto d_lb = const_local_view_device(lb);
        auto d_ub = const_local_view_device(ub);
        auto d_x = const_local_view_device(x);
        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar li = d_lb.get(i);
                       Scalar ui = d_ub.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi < li || xi > ui) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    } else if (constraints_.has_upper_bound() &&
               !constraints_.has_lower_bound()) {
      const auto &ub = *constraints_.upper_bound();

      {
        auto d_ub = const_local_view_device(ub);
        auto d_x = const_local_view_device(x);

        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar ui = d_ub.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi > ui) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    } else {
      const auto &lb = *constraints_.lower_bound();

      {
        auto d_lb = const_local_view_device(lb);
        auto d_x = const_local_view_device(x);

        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar li = d_lb.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi < li) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    }
  }

 protected:
  BoxConstraints constraints_;
  Vector help_;
};

template <class Matrix, class Vector>
class ConstrainedExtendedTestFunction
    : virtual public TestFunctionInterface<Matrix, Vector>,
      virtual public ExtendedFunction<Matrix, Vector> {
 public:
  using Scalar = typename utopia::Traits<Matrix>::Scalar;
  using SizeType = typename utopia::Traits<Vector>::SizeType;
  using BoxConstraints = utopia::BoxConstraints<Vector>;

  ~ConstrainedExtendedTestFunction() override = default;

  virtual const BoxConstraints &box_constraints() const { return constraints_; }

  virtual void set_box_constraints(const BoxConstraints &box) {
    constraints_ = box;
  }

  virtual const Vector &upper_bound() const {
    if (!constraints_.upper_bound()) {
      utopia_error("ConstrainedTestFunction::upper bound does not exist. \n");
    }

    return *constraints_.upper_bound();
  }

  virtual const Vector &lower_bound() const {
    if (!constraints_.lower_bound()) {
      utopia_error("ConstrainedTestFunction::lower bound does not exist. \n");
    }

    return *constraints_.lower_bound();
  }

  virtual bool has_lower_bound() const {
    return constraints_.has_lower_bound();
  }

  virtual bool has_upper_bound() const {
    return constraints_.has_upper_bound();
  }

  Vector initial_guess() const override = 0;

  virtual bool is_feasible(Vector &x) {
    if (!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
      return true;

    if (empty(help_) || size(help_) != size(x)) {
      help_ = 0.0 * x;
    }

    if (constraints_.has_upper_bound() && constraints_.has_lower_bound()) {
      const auto &ub = *constraints_.upper_bound();
      const auto &lb = *constraints_.lower_bound();

      {
        auto d_lb = const_local_view_device(lb);
        auto d_ub = const_local_view_device(ub);
        auto d_x = const_local_view_device(x);

        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar li = d_lb.get(i);
                       Scalar ui = d_ub.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi < li || xi > ui) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    } else if (constraints_.has_upper_bound() &&
               !constraints_.has_lower_bound()) {
      const auto &ub = *constraints_.upper_bound();

      {
        auto d_ub = const_local_view_device(ub);
        auto d_x = const_local_view_device(x);

        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar ui = d_ub.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi > ui) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    } else {
      const auto &lb = *constraints_.lower_bound();

      {
        auto d_lb = const_local_view_device(lb);
        auto d_x = const_local_view_device(x);

        auto help_view = local_view_device(help_);

        parallel_for(local_range_device(help_),
                     UTOPIA_LAMBDA(const SizeType i) {
                       Scalar li = d_lb.get(i);
                       Scalar xi = d_x.get(i);

                       help_view.set(i, (xi < li) ? 1.0 : 0.0);
                     });
      }

      return (sum(help_) > 0.0) ? false : true;
    }
  }

 protected:
  BoxConstraints constraints_;
  Vector help_;
};

}  // namespace utopia
#endif  // UTOPIA_UNCONSTRAINED_TEST_FUNCTIONS
