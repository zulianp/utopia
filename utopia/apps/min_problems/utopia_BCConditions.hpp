#ifndef UTOPIA_BC_CONDITIONS_MIN_PROBLEMS_HPP
#define UTOPIA_BC_CONDITIONS_MIN_PROBLEMS_HPP

#include "utopia_Base.hpp"
#include "utopia_RangeDevice.hpp"

// include edsl components
#include "utopia_BCSetup.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_impl.hpp"

#include <chrono>
#include <cmath>
#include <random>

namespace utopia {

template <class FunctionSpace>
class AllZeroScalar : public BCSetup<FunctionSpace>,
                      public InitialCondition<FunctionSpace> {
 public:
  using Scalar = typename FunctionSpace::Scalar;
  using Vector = typename FunctionSpace::Vector;

  AllZeroScalar(FunctionSpace &space)
      : BCSetup<FunctionSpace>(space), InitialCondition<FunctionSpace>(space) {}

  void emplace_BC() override {
    static const int Dim = FunctionSpace::Dim;
    using Point = typename FunctionSpace::Point;

    this->space_.reset_bc();

    this->space_.emplace_dirichlet_condition(
        SideSet::left(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
        0);

    this->space_.emplace_dirichlet_condition(
        SideSet::right(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
        0);

    this->space_.emplace_dirichlet_condition(
        SideSet::top(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
        0);

    this->space_.emplace_dirichlet_condition(
        SideSet::bottom(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
        0);

    if (Dim == 3) {
      this->space_.emplace_dirichlet_condition(
          SideSet::front(),
          UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; }, 0);

      this->space_.emplace_dirichlet_condition(
          SideSet::back(), UTOPIA_LAMBDA(const Point &)->Scalar { return 0.0; },
          0);
    }
  }

  void init(PetscVector &x) override {
    using Point = typename FunctionSpace::Point;
    using Dev = typename FunctionSpace::Device;
    using Mesh = typename FunctionSpace::Mesh;
    using Elem = typename FunctionSpace::Shape;
    using ElemViewScalar =
        typename utopia::FunctionSpace<Mesh, 1, Elem>::ViewDevice::Elem;
    static const int NNodes = Elem::NNodes;

    auto C = this->space_;

    auto sampler = utopia::sampler(C, UTOPIA_LAMBDA(const Point &x)->Scalar {
      Scalar f = 0.0;
      return 0.0;
    });

    {
      auto C_view = C.view_device();
      auto sampler_view = sampler.view_device();
      auto x_view = this->space_.assembly_view_device(x);

      Dev::parallel_for(this->space_.element_range(),
                        UTOPIA_LAMBDA(const SizeType &i) {
                          ElemViewScalar e;
                          C_view.elem(i, e);

                          StaticVector<Scalar, NNodes> s;
                          sampler_view.assemble(e, s);
                          C_view.set_vector(e, s, x_view);
                        });
    }
  }
};

}  // namespace utopia

#endif  // UTOPIA_BC_CONDITIONS_MIN_PROBLEMS_HPP
