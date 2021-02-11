#ifndef UTOPIA_VC_ISOTROPIC_PHASE_FIELD_HPP
#define UTOPIA_VC_ISOTROPIC_PHASE_FIELD_HPP

#include "utopia_CoefStrainView.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_ExtendedFunction.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PhaseFieldBase.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"
#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#include "utopia_IsotropicPhaseField.hpp"

#include "utopia_AppBase.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

// #define ENABLE_ASTRUM_CONDITIONS

#ifdef USE_SIMD_ASSEMBLY
#define USE_SIMD_PHASE_FIELD
#endif

#ifdef USE_SIMD_PHASE_FIELD
namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class VcIsotropicPhaseFieldForBrittleFractures final : public PhaseFieldFracBase<FunctionSpace, Dim> {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;
        using Parameters = utopia::PFFracParameters<FunctionSpace>;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        // using Quadrature = utopia::Quadrature<Shape, 2*(Shape::Order -1)>;
        // using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        using SIMDType = Vc::Vector<Scalar>;
        using Quadrature = simd_v2::Quadrature<Scalar, Dim>;
        using CGradValue = typename simd_v2::FETraits<CElem, Scalar>::GradValue;
        using UGradValue = typename simd_v2::FETraits<UElem, Scalar>::GradValue;
        // using GradValue = typename simd_v2::FETraits<Elem, SIMDType>::GradValue;
        // using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        VcIsotropicPhaseFieldForBrittleFractures(FunctionSpace &space)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space) {}

        VcIsotropicPhaseFieldForBrittleFractures(FunctionSpace &space, const Parameters &params)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space, params) {}

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("VcIsotropicPhaseFieldForBrittleFractures::value");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace> >(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;
            simd_v2::QuadratureDB<Shape>::get(2 * Shape::Order, q);

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_fun.coefficient(), q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                const SizeType n_qp = q.n_points();
                Vc::vector<SIMDType> c(n_qp, 0.0);
                Vc::vector<SIMDType> c_old(n_qp, 0.0);
                Vc::vector<SIMDType> p(n_qp, 0.0);
                Vc::vector<CGradValue> c_grad_el(n_qp);
                Vc::vector<UGradValue> el_strain(n_qp);

                Device::parallel_reduce(
                    this->space_.element_range(),
                    [&](const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);
                        p_view.get(c_e, p);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        strain_view.get(u_e, el_strain);
                        c_grad_view.get(c_e, c_grad_el);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;
                        for (SizeType qp = 0; qp < n_qp; ++qp) {
                            auto dx_qp = dx(qp);
                            auto tr = trace(el_strain[qp]);
                            if (this->params_.use_pressure) {
                                el_energy +=
                                    simd_v2::integrate(PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(
                                                           this->params_, c[qp]) *
                                                       p[qp] * tr * dx_qp);
                            }

                            el_energy += simd_v2::integrate(
                                energy(this->params_, c[qp], c_grad_el[qp], tr, el_strain[qp]) * dx_qp);

                            if (this->params_.use_penalty_irreversibility) {
                                auto c_cold = c[qp] - c_old[qp];
                                auto c_mask = c_cold < 0.0;
                                auto c_cold_bracket = c_cold;
                                // auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                c_cold_bracket.setZeroInverted(c_mask);

                                el_energy += simd_v2::integrate(this->params_.penalty_param / 2.0 * c_cold_bracket *
                                                                c_cold_bracket * dx_qp);
                            }
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            if (!empty(this->force_field_)) {
                // MAYBE -= dot(x_const, this->force_field_);
                val += dot(x_const, this->force_field_);
            }

            // this->add_pf_constraints(x_const);

            UTOPIA_TRACE_REGION_END("VcIsotropicPhaseFieldForBrittleFractures::value");
            return true;
        }

        bool elastic_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("VcIsotropicPhaseFieldForBrittleFractures::elastic_energy");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace> >(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;
            simd_v2::QuadratureDB<Shape>::get(2 * Shape::Order, q);

            ////////////////////////////////////////////////////////////////////////////

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                const SizeType n_qp = q.n_points();
                Vc::vector<SIMDType> c(n_qp, 0.0);
                Vc::vector<SIMDType> c_old(n_qp, 0.0);
                Vc::vector<UGradValue> el_strain(n_qp);

                Device::parallel_reduce(
                    this->space_.element_range(),
                    [&](const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        strain_view.get(u_e, el_strain);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;

                        for (SizeType qp = 0; qp < n_qp; ++qp) {
                            auto tr = trace(el_strain[qp]);

                            el_energy +=
                                simd_v2::integrate(elastic_energy(this->params_, c[qp], tr, el_strain[qp]) * dx(qp));
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            UTOPIA_TRACE_REGION_END("VcIsotropicPhaseFieldForBrittleFractures::elastic_energy");
            return true;
        }

        bool fracture_energy(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("VcIsotropicPhaseFieldForBrittleFractures::fracture_energy");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace> >(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;
            simd_v2::QuadratureDB<Shape>::get(2 * Shape::Order, q);

            ////////////////////////////////////////////////////////////////////////////

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();

                auto differential_view = differential.view_device();

                const SizeType n_qp = q.n_points();
                Vc::vector<SIMDType> c(n_qp, 0.0);
                Vc::vector<SIMDType> c_old(n_qp, 0.0);
                Vc::vector<CGradValue> c_grad_el(n_qp);

                Device::parallel_reduce(
                    this->space_.element_range(),
                    [&](const SizeType &i) {
                        CElem c_e;
                        C_view.elem(i, c_e);

                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);
                        c_grad_view.get(c_e, c_grad_el);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;

                        for (SizeType qp = 0; qp < n_qp; ++qp) {
                            el_energy += simd_v2::integrate(PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
                                                                this->params_, c[qp], c_grad_el[qp]) *
                                                            dx(qp));
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            UTOPIA_TRACE_REGION_END("VcIsotropicPhaseFieldForBrittleFractures::fracture_energy");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("VcIsotropicPhaseFieldForBrittleFractures::gradient");

            if (empty(g)) {
                this->space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace> >(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;
            simd_v2::QuadratureDB<Shape>::get(2 * Shape::Order, q);

            ////////////////////////////////////////////////////////////////////////////

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto press_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            CoefStrain<USpace, Quadrature> strain(u_coeff, q);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_old_view = c_old.view_device();
                auto p_view = press_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = this->space_.assembly_view_device(g);
                auto ref_strain_u_view = ref_strain_u.view_device();

                const SizeType n_qp = q.n_points();
                Vc::vector<SIMDType> c(n_qp, 0.0);
                Vc::vector<SIMDType> c_old(n_qp, 0.0);
                Vc::vector<SIMDType> p(n_qp);
                Vc::vector<CGradValue> c_grad_el(n_qp);
                Vc::vector<UGradValue> el_strain(n_qp);

                simd_v2::Matrix<Scalar, Dim, Dim> stress;
                StaticVector<Scalar, U_NDofs> u_el_vec;
                StaticVector<Scalar, C_NDofs> c_el_vec;

                Device::parallel_for(this->space_.element_range(), [&](const SizeType &i) {
                    u_el_vec.set(0.0);
                    c_el_vec.set(0.0);

                    ////////////////////////////////////////////

                    UElem u_e;
                    U_view.elem(i, u_e);
                    strain_view.get(u_e, el_strain);
                    // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                    auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                    ////////////////////////////////////////////

                    CElem c_e;
                    C_view.elem(i, c_e);

                    c_view.get(c_e, c);

                    c_old_view.get(c_e, c_old);

                    p_view.get(c_e, p);

                    c_grad_view.get(c_e, c_grad_el);
                    auto dx = differential_view.make(c_e);
                    auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                    auto c_shape_fun_el = c_shape_view.make(c_e);

                    ////////////////////////////////////////////

                    for (SizeType qp = 0; qp < n_qp; ++qp) {
                        const auto dx_qp = dx(qp);
                        auto tr_strain_u = trace(el_strain[qp]);

                        compute_stress(this->params_, tr_strain_u, el_strain[qp], stress);
                        auto gc_qp =
                            PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(this->params_, c[qp]);
                        stress = (gc_qp * (1.0 - this->params_.regularization) + this->params_.regularization) * stress;

                        for (SizeType j = 0; j < U_NDofs; ++j) {
                            auto &&strain_test = u_strain_shape_el(j, qp);
                            u_el_vec(j) += simd_v2::integrate(inner(stress, strain_test) * dx_qp);

                            if (this->params_.use_pressure) {
                                u_el_vec(j) += simd_v2::integrate(gc_qp * p[qp] * trace(strain_test) * dx_qp);
                            }
                        }

                        const auto elast = grad_elastic_energy_wrt_c(this->params_, c[qp], tr_strain_u, el_strain[qp]);

                        for (SizeType j = 0; j < C_NDofs; ++j) {
                            const auto shape_test = c_shape_fun_el(j, qp);
                            const auto frac = PhaseFieldFracBase<FunctionSpace, Dim>::grad_fracture_energy_wrt_c(
                                this->params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                            if (this->params_.use_pressure) {
                                const auto der_c_pres =
                                    PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(this->params_,
                                                                                                        c[qp]) *
                                    p[qp] * tr_strain_u * shape_test;
                                c_el_vec(j) += simd_v2::integrate(der_c_pres * dx_qp);
                            }

                            c_el_vec(j) += simd_v2::integrate((elast * shape_test + frac) * dx_qp);

                            if (this->params_.use_penalty_irreversibility) {
                                auto c_cold = c[qp] - c_old[qp];
                                auto c_mask = c_cold < 0.0;
                                auto c_cold_bracket = c_cold;
                                // auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                c_cold_bracket.setZeroInverted(c_mask);

                                c_el_vec(j) += simd_v2::integrate(this->params_.penalty_param * c_cold_bracket *
                                                                  shape_test * dx_qp);
                            }
                        }
                    }

                    U_view.add_vector(u_e, u_el_vec, g_view);
                    C_view.add_vector(c_e, c_el_vec, g_view);
                });
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_grad(*this, x_const, g);
            }

            if (!empty(this->force_field_)) {
                // MAYBE g -= this->force_field_;
                g += this->force_field_;
            }

            this->space_.apply_zero_constraints(g);

            // fully broken case is treated as Dirichlet BC
            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(g, x_const);
            }

            UTOPIA_TRACE_REGION_END("VcIsotropicPhaseFieldForBrittleFractures::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("VcIsotropicPhaseFieldForBrittleFractures::hessian");

            if (empty(H)) {
                // if(use_dense_hessian_) {
                //     H = local_zeros({this->space_.n_dofs(), this->space_.n_dofs()}); //FIXME
                // } else {
                this->space_.create_matrix(H);
                // }
            } else {
                H *= 0.0;
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            // auto &x     = const_cast<Vector &>(x_const);t
            // auto &press = const_cast<Vector &>(pressure_field_);

            ////////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace> >(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *this->local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_pressure_field_);

            // update c_old
            this->space_.global_to_local(this->x_old_, *this->local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace> >(C, this->local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;
            simd_v2::QuadratureDB<Shape>::get(2 * Shape::Order, q);

            ////////////////////////////////////////////////////////////////////////////

            auto c_val = c_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            // value based
            CoefStrain<USpace, Quadrature> strain(u_coeff, q);

            // reference based
            ShapeStress<USpace, Quadrature> p_stress(U, q, this->params_.mu, this->params_.lambda);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = this->space_.view_device();

                auto c_view = c_val.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();
                auto p_stress_view = p_stress.view_device();

                auto H_view = this->space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                const SizeType n_qp = q.n_points();
                Vc::vector<SIMDType> c(n_qp, 0.0);
                Vc::vector<SIMDType> p(n_qp, 0.0);
                Vc::vector<UGradValue> el_strain(n_qp);

                StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;

                StaticMatrix<Scalar, U_NDofs, U_NDofs> u_mat;
                StaticMatrix<Scalar, C_NDofs, C_NDofs> c_mat;
                StaticMatrix<Scalar, C_NDofs, U_NDofs> cu_mat;

                simd_v2::Matrix<Scalar, Dim, Dim> stress;

                Device::parallel_for(this->space_.element_range(), [&](const SizeType &i) {
                    MixedElem e;
                    space_view.elem(i, e);
                    // el_mat.set(0.0);

                    u_mat.set(0.0);
                    c_mat.set(0.0);
                    cu_mat.set(0.0);

                    ////////////////////////////////////////////
                    UElem u_e;
                    U_view.elem(i, u_e);
                    strain_view.get(u_e, el_strain);
                    // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                    auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                    ////////////////////////////////////////////
                    CElem c_e;
                    C_view.elem(i, c_e);
                    c_view.get(c_e, c);
                    p_view.get(c_e, p);

                    auto dx = differential_view.make(c_e);
                    auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                    auto c_shape_fun_el = c_shape_view.make(c_e);

                    ////////////////////////////////////////////
                    for (SizeType qp = 0; qp < n_qp; ++qp) {
                        const auto tr_strain_u = trace(el_strain[qp]);
                        const auto dx_qp = dx(qp);

                        const auto eep = elastic_energy(this->params_, c[qp], tr_strain_u, el_strain[qp]);

                        // pragma GCCunroll(C_NDofs)
                        for (SizeType l = 0; l < C_NDofs; ++l) {
                            const auto c_shape_l = c_shape_fun_el(l, qp);
                            auto &&c_grad_l = c_grad_shape_el(l, qp);

                            // SYMMETRIC VERSION
                            for (SizeType j = l; j < C_NDofs; ++j) {
                                const auto c_shape_j_l_prod = c_shape_fun_el(j, qp) * c_shape_l;

                                auto val =
                                    bilinear_cc(
                                        this->params_, c[qp], eep, c_shape_j_l_prod, c_grad_shape_el(j, qp), c_grad_l) *
                                    dx_qp;

                                if (this->params_.use_pressure) {
                                    val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv2(
                                               this->params_, c[qp]) *
                                           p[qp] * tr_strain_u * c_shape_j_l_prod * dx_qp;
                                }

                                if (this->params_.use_penalty_irreversibility) {
                                    val += this->params_.penalty_param * c_shape_j_l_prod * dx_qp;
                                }

                                auto val_integr = simd_v2::integrate((l == j) ? (0.5 * val) : val);

                                // el_mat(l, j) += val_integr;
                                // el_mat(j, l) += val_integr;

                                c_mat(l, j) += val_integr;
                                c_mat(j, l) += val_integr;
                            }
                        }

                        for (SizeType l = 0; l < U_NDofs; ++l) {
                            auto &&u_strain_shape_l = u_strain_shape_el(l, qp);

                            for (SizeType j = l; j < U_NDofs; ++j) {
                                auto val = PhaseFieldFracBase<FunctionSpace, Dim>::bilinear_uu(
                                               this->params_, c[qp], p_stress_view(j, qp), u_strain_shape_l) *
                                           dx_qp;

                                auto val_integr = simd_v2::integrate((l == j) ? (0.5 * val) : val);
                                // el_mat(C_NDofs + l, C_NDofs + j) += val_integr;
                                // el_mat(C_NDofs + j, C_NDofs + l) += val_integr;

                                u_mat(l, j) += val_integr;
                                u_mat(j, l) += val_integr;
                            }
                        }

                        //////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_ASTRUM_CONDITIONS
                        if (!this->params_.turn_off_uc_coupling || !this->params_.turn_off_cu_coupling) {
#endif  // ENABLE_ASTRUM_CONDITIONS
                            compute_stress(this->params_, tr_strain_u, el_strain[qp], stress);
                            for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                                // CHANGE (pre-compute/store shape fun)
                                const auto c_shape_i = c_shape_fun_el(c_i, qp);

                                for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                                    auto &&strain_shape = u_strain_shape_el(u_i, qp);

                                    auto val =
                                        bilinear_uc(this->params_, c[qp], stress, strain_shape, c_shape_i) * dx_qp;

                                    if (this->params_.use_pressure) {
                                        const auto tr_strain_shape = trace(strain_shape);
                                        val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(
                                                   this->params_, c[qp]) *
                                               p[qp] * tr_strain_shape * c_shape_i * dx_qp;
                                    }

                                    // not symetric, but more numerically stable
                                    if (this->params_.turn_off_cu_coupling == false) {
                                        // el_mat(c_i, C_NDofs + u_i) += simd_v2::integrate(val);
                                        cu_mat(c_i, u_i) += simd_v2::integrate(val);
                                    }

                                    if (this->params_.turn_off_uc_coupling == false) {
                                        // el_mat(C_NDofs + u_i, c_i) += simd_v2::integrate(val);
                                        cu_mat(C_NDofs + u_i, c_i) += simd_v2::integrate(val);
                                    }
                                }
                            }
#ifdef ENABLE_ASTRUM_CONDITIONS
                        }
#endif  // ENABLE_ASTRUM_CONDITIONS
                    }

                    el_mat.set_matrix(0, 0, c_mat);
                    el_mat.set_matrix(C_NDofs, C_NDofs, u_mat);
                    el_mat.set_matrix(C_NDofs, 0, transpose(cu_mat));
                    el_mat.set_matrix(0, C_NDofs, cu_mat);

                    space_view.add_matrix(e, el_mat, H_view);
                });
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_hessian(*this, x_const, H);
            }

            this->space_.apply_constraints(H);

            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(H, x_const);
            }

            UTOPIA_TRACE_REGION_END("VcIsotropicPhaseFieldForBrittleFractures::hessian");
            return true;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        template <typename PhaseFieldValue, typename ElasticEnergy, typename ShapeProd, class GradShape>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_cc(const Parameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const ElasticEnergy &elastic_energy,
                                                                  const ShapeProd &shape_prod,
                                                                  const GradShape &grad_trial,
                                                                  const GradShape &grad_test) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::diffusion_c(params, grad_trial, grad_test) +
                   reaction_c(params, shape_prod) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy, shape_prod);
        }

        // (sigma+(phi_u), epsilon(u)) * g'_c * phi_c
        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cu(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return ((1.0 - params.regularization) *
                    PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value)) *
                   c_trial_fun * inner(stress_p, full_strain);
        }

        template <typename PhaseFieldValue, class Stress, class FullStrain, typename CTrialFun>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uc(const Parameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const Stress &stress,
                                                                  const FullStrain &full_strain,
                                                                  const CTrialFun &c_trial_fun) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value) *
                   c_trial_fun * inner(stress, full_strain);
        }

        template <typename ShapeProd>
        UTOPIA_INLINE_FUNCTION static ShapeProd reaction_c(const Parameters &params, const ShapeProd &shape_prod) {
            return (params.fracture_toughness / params.length_scale) * shape_prod;
        }

        template <typename PhaseFieldValue, typename ElasticEnergy, typename ShapeProd>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue elastic_deriv_cc(const Parameters &params,
                                                                       const PhaseFieldValue &phase_field_value,
                                                                       const ElasticEnergy &elastic_energy,
                                                                       const ShapeProd &shape_prod) {
            const PhaseFieldValue dcc =
                (1.0 - params.regularization) *
                PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * shape_prod * elastic_energy;
        }

        template <typename PhaseFiledValue, typename TraceT, class Strain>
        UTOPIA_INLINE_FUNCTION static PhaseFiledValue grad_elastic_energy_wrt_c(
            const Parameters &params,
            const PhaseFiledValue &phase_field_value,
            const TraceT &trace,
            const Strain &strain) {
            return (PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value) *
                    (1.0 - params.regularization)) *
                   strain_energy(params, trace, strain);
        }

        template <typename TraceT, class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const Parameters &params,
                                                          const TraceT &tr,
                                                          const Strain &strain,
                                                          Stress &stress) {
            stress = (2.0 * params.mu * strain) + (params.lambda * tr * (device::identity<Scalar>()));
        }

        template <class Grad, class PhaseFieldValue, class Strain>
        UTOPIA_INLINE_FUNCTION static SIMDType energy(const Parameters &params,
                                                      // c
                                                      const PhaseFieldValue &phase_field_value,
                                                      const Grad &phase_field_grad,
                                                      // u
                                                      const SIMDType &trace,
                                                      const Strain &strain) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
                       params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, trace, strain);
        }

        template <typename TraceT, class Strain>
        UTOPIA_INLINE_FUNCTION static TraceT strain_energy(const Parameters &params,
                                                           const TraceT trace,
                                                           const Strain &strain) {
            return 0.5 * params.lambda * trace * trace + params.mu * inner(strain, strain);
        }

        template <typename PhaseFieldValue, typename TraceT, class Strain>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue elastic_energy(const Parameters &params,
                                                                     const PhaseFieldValue &phase_field_value,
                                                                     const TraceT &trace,
                                                                     const Strain &strain) {
            return (PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(params, phase_field_value) *
                        (1.0 - params.regularization) +
                    params.regularization) *
                   strain_energy(params, trace, strain);
        }
    };

}  // namespace utopia

#else
namespace utopia {
    template <class FunctionSpace>
    using VcIsotropicPhaseFieldForBrittleFractures = IsotropicPhaseFieldForBrittleFractures<FunctionSpace>;
}
#endif  // USE_SIMD_PHASE_FIELD

// clean-up macros
#undef ENABLE_ASTRUM_CONDITIONS
#undef UNROLL_FACTOR
#undef U_MIN
#endif
