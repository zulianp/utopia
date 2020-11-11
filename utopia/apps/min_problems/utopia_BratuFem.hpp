#ifndef UTOPIA_QP_POISSON_HPP
#define UTOPIA_QP_POISSON_HPP

#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"

namespace utopia {

template <class FunctionSpace, int Dim = FunctionSpace::Dim>
// class BratuFem final : public ExtendedFunction<typename
// FunctionSpace::Matrix, typename FunctionSpace::Vector>
class BratuFem final
    : virtual public UnconstrainedExtendedTestFunction<
          typename FunctionSpace::Matrix, typename FunctionSpace::Vector>,
      virtual public ConstrainedExtendedTestFunction<
          typename FunctionSpace::Matrix, typename FunctionSpace::Vector> {
 public:
  using Scalar = typename FunctionSpace::Scalar;
  using SizeType = typename FunctionSpace::SizeType;
  using Vector = typename FunctionSpace::Vector;
  using Matrix = typename FunctionSpace::Matrix;
  using Device = typename FunctionSpace::Device;

  using CSpace = typename FunctionSpace::template Subspace<1>;
  using CElem = typename CSpace::ViewDevice::Elem;

  // FIXME
  using Shape = typename FunctionSpace::Shape;
  using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

  static const int C_NDofs = CSpace::NDofs;
  static const int NQuadPoints = Quadrature::NPoints;

  class Parameters : public Configurable {
   public:
    void read(Input &in) override {
      in.get("lambda", lambda);
      in.get("rhs_type", rhs_type);
    }

    Parameters() : lambda(3.0), rhs_type(0) {}

    Scalar lambda;
    SizeType rhs_type;
  };

  void read(Input &in) override {
    params_.read(in);

    this->init_forcing_function(this->rhs_);
  }

  BratuFem(FunctionSpace &space) : space_(space) {
    // needed for ML setup
    space_.create_vector(this->_x_eq_values);
    space_.create_vector(this->_eq_constrains_flg);
    space_.create_vector(this->x_exact_);
    x_exact_.set(0.0);

    this->local_x_ = std::make_shared<Vector>();
    space_.create_local_vector(*this->local_x_);

    // TODO:: add more options
    this->init_forcing_function(this->rhs_);
    this->init_constraints();
  }

  inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
    space_.create_matrix(H);
    return true;
  }

  inline bool update(const Vector & /*x*/) override {
    // x_coeff_.update(x);
    return true;
  }

  bool value(const Vector &x_const, Scalar &val) const override {
    UTOPIA_TRACE_REGION_BEGIN("QPPoisson::value");

    CSpace C = space_.subspace(0);

    // update local vector x
    space_.global_to_local(x_const, *local_x_);
    auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

    FEFunction<CSpace> c_fun(c_coeff);
    ////////////////////////////////////////////////////////////////////////////
    Quadrature q;

    auto c_val = c_fun.value(q);
    auto c_grad = c_fun.gradient(q);
    auto differential = C.differential(q);
    auto c_shape = C.shape(q);

    val = 0.0;

    {
      auto C_view = C.view_device();
      auto c_view = c_val.view_device();
      auto c_grad_view = c_grad.view_device();
      auto differential_view = differential.view_device();

      Device::parallel_reduce(
          space_.element_range(),
          UTOPIA_LAMBDA(const SizeType &i) {

            CElem c_e;
            C_view.elem(i, c_e);

            StaticVector<Scalar, NQuadPoints> c;
            c_view.get(c_e, c);

            auto c_grad_el = c_grad_view.make(c_e);
            auto dx = differential_view.make(c_e);

            Scalar el_energy = 0.0;
            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
              el_energy += 0.5 * inner(c_grad_el[qp], c_grad_el[qp]) * dx(qp);

              el_energy -= params_.lambda * device::exp(c[qp]) * dx(qp);
            }

            assert(el_energy == el_energy);
            return el_energy;
          },
          val);
    }

    val = x_const.comm().sum(val);

    // add contributions from RHS
    val = val + dot(x_const, this->rhs_);

    UTOPIA_TRACE_REGION_END("QPPoisson::value");
    return true;
  }

  bool gradient(const Vector &x_const, Vector &g) const override {
    UTOPIA_TRACE_REGION_BEGIN("QPPoisson::gradient");

    if (empty(g)) {
      space_.create_vector(g);
    } else {
      g.set(0.0);
    }

    CSpace C = space_;

    ///////////////////////////////////////////////////////////////////////////

    // update local vector x
    space_.global_to_local(x_const, *local_x_);
    auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

    FEFunction<CSpace> c_fun(c_coeff);

    Quadrature q;
    auto c_val = c_fun.value(q);
    auto c_grad = c_fun.gradient(q);
    auto differential = C.differential(q);

    auto c_shape = C.shape(q);
    auto c_grad_shape = C.shape_grad(q);

    {
      auto C_view = C.view_device();
      auto c_view = c_val.view_device();
      auto c_grad_view = c_grad.view_device();

      auto differential_view = differential.view_device();

      // auto v_grad_shape_view = v_grad_shape.view_device();
      auto c_shape_view = c_shape.view_device();
      auto c_grad_shape_view = c_grad_shape.view_device();

      auto g_view = space_.assembly_view_device(g);

      Device::parallel_for(space_.element_range(), UTOPIA_LAMBDA(
                                                       const SizeType &i) {
        StaticVector<Scalar, C_NDofs> c_el_vec;
        c_el_vec.set(0.0);
        ////////////////////////////////////////////

        CElem c_e;
        C_view.elem(i, c_e);
        StaticVector<Scalar, NQuadPoints> c;
        c_view.get(c_e, c);

        auto c_grad_el = c_grad_view.make(c_e);
        auto dx = differential_view.make(c_e);
        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
        auto c_shape_fun_el = c_shape_view.make(c_e);

        ////////////////////////////////////////////
        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
          for (SizeType j = 0; j < C_NDofs; ++j) {
            c_el_vec(j) +=
                inner(c_grad_el[qp], c_grad_shape_el(j, qp)) * dx(qp);

            const Scalar shape_test = c_shape_fun_el(j, qp);
            c_el_vec(j) -=
                inner(params_.lambda * device::exp(c[qp]), shape_test) * dx(qp);
          }
        }

        C_view.add_vector(c_e, c_el_vec, g_view);
      });
    }

    g = g + this->rhs_;
    space_.apply_zero_constraints(g);

    UTOPIA_TRACE_REGION_END("QPPoisson::gradient");
    return true;
  }

  bool hessian(const Vector &x_const, Matrix &H) const override {
    UTOPIA_TRACE_REGION_BEGIN("QPPoisson::hessian");

    if (empty(H)) {
      // if(use_dense_hessian_) {
      //     H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
      // } else {
      space_.create_matrix(H);
      // }
    } else {
      H *= 0.0;
    }

    CSpace C = space_;

    ////////////////////////////////////////////////////////////////////////////
    // update local vector x
    space_.global_to_local(x_const, *local_x_);
    auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

    FEFunction<CSpace> c_fun(c_coeff);

    ////////////////////////////////////////////////////////////////////////////

    Quadrature q;

    auto c_val = c_fun.value(q);
    auto c_grad = c_fun.gradient(q);
    auto differential = C.differential(q);

    auto c_shape = C.shape(q);
    auto c_grad_shape = C.shape_grad(q);

    {
      auto C_view = C.view_device();

      auto space_view = space_.view_device();

      auto c_view = c_val.view_device();
      auto c_grad_view = c_grad.view_device();
      auto differential_view = differential.view_device();

      auto c_shape_view = c_shape.view_device();
      auto c_grad_shape_view = c_grad_shape.view_device();

      auto H_view = space_.assembly_view_device(H);

      Device::parallel_for(
          space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
            StaticMatrix<Scalar, C_NDofs, C_NDofs> el_mat;
            el_mat.set(0.0);

            ////////////////////////////////////////////
            CElem c_e;
            C_view.elem(i, c_e);
            StaticVector<Scalar, NQuadPoints> c;
            c_view.get(c_e, c);

            auto dx = differential_view.make(c_e);
            auto c_grad_shape_el = c_grad_shape_view.make(c_e);
            auto c_shape_fun_el = c_shape_view.make(c_e);

            ////////////////////////////////////////////
            for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
              for (SizeType l = 0; l < C_NDofs; ++l) {
                auto &&c_grad_l = c_grad_shape_el(l, qp);
                const Scalar c_shape_l = c_shape_fun_el(l, qp);

                Scalar val1 = 0;
                val1 -= inner(params_.lambda * device::exp(c[qp]), c_shape_l) *
                        dx(qp);
                el_mat(l, l) += val1;

                for (SizeType j = l; j < C_NDofs; ++j) {
                  Scalar val = inner(c_grad_shape_el(j, qp), c_grad_l) * dx(qp);
                  // val -= inner(params_.lambda * device::exp(c[qp]),
                  //              c_shape_fun_el(l, qp)) *
                  //        dx(qp);

                  val = (l == j) ? (0.5 * val) : val;
                  el_mat(l, j) += val;
                  el_mat(j, l) += val;
                }
              }
            }

            space_view.add_matrix(c_e, el_mat, H_view);
          });
    }

    space_.apply_constraints(H);

    UTOPIA_TRACE_REGION_END("QPPoisson::hessian");
    return true;
  }

 private:
  void init_forcing_function(Vector &rhs) {
    using Elem = typename FunctionSpace::Elem;
    static const int NNodes = Elem::NNodes;
    using ElementVector = utopia::StaticVector<Scalar, NNodes>;
    using Mesh = typename FunctionSpace::Mesh;
    using Point = typename Mesh::Point;

    if (params_.rhs_type == 0) {
      space_.create_vector(rhs);
      rhs = 0.0 * rhs;
      return;
    }

    auto forcing_function = UTOPIA_LAMBDA(const Point &x)->Scalar {
      Scalar pi = 3.1415926;

      if (params_.rhs_type == 1) {
        return 1.0;
      } else if (params_.rhs_type == 2) {
        return -4.0 * pi * pi * device::sin(2.0 * pi * x[0]) *
               device::sin(2.0 * pi * x[1]);
      } else if (params_.rhs_type == 3) {
        return -2.0 * pi * pi * device::sin(pi * x[0]) * device::sin(pi * x[1]);
      } else {
        std::cerr << "---------- forcing function choice wrong ------- \n";
        return 0.0;
      }
    };

    space_.create_vector(rhs);
    Quadrature q;

    Projection<FunctionSpace, Quadrature, decltype(forcing_function)> proj(
        space_, q, forcing_function);

    auto proj_view = proj.view_device();

    {
      auto space_view = space_.view_device();
      auto rhs_view = space_.assembly_view_device(rhs);

      Device::parallel_for(space_.element_range(),
                           UTOPIA_LAMBDA(const SizeType &i) {
                             Elem e;
                             space_view.elem(i, e);

                             ElementVector el_vec;
                             el_vec.set(0.0);

                             proj_view.assemble(i, e, el_vec);
                             space_view.add_vector(e, el_vec, rhs_view);
                           });
    }
  }

  void init_constraints() {
    // ToDO:: fix based on example
    Vector lb;
    space_.create_vector(lb);
    lb.set(-1.0);
    this->constraints_ =
        make_lower_bound_constraints(std::make_shared<Vector>(lb));
  }

 public:
  Vector initial_guess() const override {
    Vector solution;
    space_.create_vector(solution);
    rename("X", solution);

    // TODO:: add initial condition
    solution.set(-1.0);
    space_.apply_constraints(solution);

    return solution;
  }

  // not known...
  const Vector &exact_sol() const {
    std::cout << "BratuFEM:: exact Solution not known, terminate... \n";
    return x_exact_;
  }

  Scalar min_function_value() const {
    std::cout << "BratuFEM:: min_function_value not known, terminate... \n";
    return 0;
  }

  virtual std::string name() const { return "BratuFem"; }

  virtual SizeType dim() const { return size(rhs_); }

  virtual bool exact_sol_known() const { return false; }

  virtual bool parallel() const { return true; }

 private:
  FunctionSpace &space_;
  Parameters params_;

  Vector rhs_;
  Vector x_exact_;

  std::shared_ptr<Vector> local_x_;
};

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
