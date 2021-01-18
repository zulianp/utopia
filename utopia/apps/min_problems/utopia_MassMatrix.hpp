#ifndef UTOPIA_MASS_MATRIX_ONE_VAR_HPP
#define UTOPIA_MASS_MATRIX_ONE_VAR_HPP

#include "utopia_DiffController.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_Views.hpp"

namespace utopia {

template <class FunctionSpace, int Dim = FunctionSpace::Dim>
class Reaction final : public Configurable

{
 public:
  using Scalar = typename FunctionSpace::Scalar;
  using SizeType = typename FunctionSpace::SizeType;
  using Vector = typename FunctionSpace::Vector;
  using Matrix = typename FunctionSpace::Matrix;
  using Device = typename FunctionSpace::Device;

  using CSpace = typename FunctionSpace::template Subspace<1>;
  using CElem = typename CSpace::ViewDevice::Elem;
  using MixedElem = typename FunctionSpace::ViewDevice::Elem;

  // FIXME
  using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;

  static const int C_NDofs = CSpace::NDofs;
  static const int NQuadPoints = Quadrature::NPoints;

  void read(Input &) override {}

  Reaction(FunctionSpace &space) : space_(space) {}

  bool mass_matrix(Matrix &H) const {
    if (empty(H)) {
      space_.create_matrix(H);
    } else {
      H *= 0.0;
    }

    CSpace C = space_.subspace(0);

    Quadrature q;

    auto differential = C.differential(q);
    auto c_shape = C.shape(q);

    {
      auto C_view = C.view_device();

      auto space_view = space_.view_device();
      auto differential_view = differential.view_device();
      auto c_shape_view = c_shape.view_device();

      auto H_view = space_.assembly_view_device(H);

      Device::parallel_for(
          space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
            StaticMatrix<Scalar, C_NDofs, C_NDofs> el_mat;
            el_mat.set(0.0);

            MixedElem e;
            space_view.elem(i, e);

            CElem c_e;
            C_view.elem(i, c_e);

            auto c_shape_fun_el = c_shape_view.make(c_e);
            auto dx = differential_view.make(c_e);

            ////////////////////////////////////////////
            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
              for (SizeType l = 0; l < C_NDofs; ++l) {
                const Scalar c_shape_l = c_shape_fun_el(l, qp);

                for (SizeType j = 0; j < C_NDofs; ++j) {
                  auto val = c_shape_fun_el(j, qp) * c_shape_l * dx(qp);
                  // c-component
                  el_mat(l, j) += val;
                }
              }
            }

            space_view.add_matrix(e, el_mat, H_view);
          });
    }

    return true;
  }

 private:
  FunctionSpace &space_;
};

}  // namespace utopia
#endif  // UTOPIA_MASS_MATRIX_ONE_VAR_HPP