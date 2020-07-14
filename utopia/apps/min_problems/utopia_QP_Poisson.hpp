#ifndef UTOPIA_QP_POISSON_HPP
#define UTOPIA_QP_POISSON_HPP

#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"


#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class QPPoisson final : public ExtendedFunction<typename FunctionSpace::Matrix, typename FunctionSpace::Vector> 
    {
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

                // in.get("mu", mu);
                // in.get("lambda", lambda);

            }

            Parameters()
                // : 
                  // a(1.0),
                  // b(1.0)
            {

            }

            // Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda;
            // Scalar regularization, pressure, penalty_param, crack_set_tol;
            // bool use_penalty_irreversibility{false}, use_crack_set_irreversibiblity{false}, use_pressure{false};
        };

        void read(Input &in) override {
            params_.read(in);

            // in.get("use_dense_hessian", use_dense_hessian_);
            // in.get("check_derivatives", check_derivatives_);
            // in.get("diff_controller", diff_ctrl_);
            
            // init_force_field(in);
        }

        QPPoisson(FunctionSpace &space) : space_(space)
        {
            // needed for ML setup
            space_.create_vector(this->_x_eq_values);
            space_.create_vector(this->_eq_constrains_flg);

            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);
        }

        QPPoisson(FunctionSpace &space, const Parameters &params) :  space_(space), params_(params)
        {
            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_pressure_field_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_pressure_field_);

            this->local_c_old_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_c_old_);
        }


        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override 
        {
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
                            el_energy += inner(c_grad_el[qp], c_grad_el[qp]) * dx(qp);
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);


            UTOPIA_TRACE_REGION_END("QPPoisson::value");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("QPPoisson::gradient");

            UTOPIA_TRACE_REGION_END("QPPoisson::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("QPPoisson::hessian");


            UTOPIA_TRACE_REGION_END("QPPoisson::hessian");
            return true;
        }


    private:
        FunctionSpace &space_;
        Parameters params_;

        std::shared_ptr<Vector> local_x_;
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
