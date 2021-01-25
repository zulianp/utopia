#ifndef UTOPIA_PHASE_FIELD_VOL_DEV_SPLIT_HPP
#define UTOPIA_PHASE_FIELD_VOL_DEV_SPLIT_HPP

#include "utopia_CoefStrainView.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PhaseFieldBase.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    // optimize 4th order tensor com
    // add pressure
    // add kappa-reg.

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PhaseFieldVolDevSplit final : public PhaseFieldFracBase<FunctionSpace, Dim>

    {
    public:
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;
        using PFFracParameters = utopia::PFFracParameters<FunctionSpace>;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;

        // FIXME
        using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;
        // using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 0>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        PhaseFieldVolDevSplit(FunctionSpace &space) : PhaseFieldFracBase<FunctionSpace, Dim>(space) {
            this->params_.fill_in_isotropic_elast_tensor();
        }

        PhaseFieldVolDevSplit(FunctionSpace &space, const PFFracParameters &params)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space, params) {
            this->params_.fill_in_isotropic_elast_tensor();
        }

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldVolDevSplit::value(...)");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            ///////////////////////////////////////////////////////////////////////////
            // update local vector x
            this->space_.global_to_local(x_const, *(this->local_x_));
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *(this->local_pressure_field_));
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto p_val = press_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;
            PrincipalStrains<USpace, Quadrature> strain(u_fun.coefficient(), q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();
                auto p_view = p_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> strain_n;
                        StaticMatrix<Scalar, Dim, Dim> strain_p;

                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> p;
                        c_view.get(c_e, c);
                        p_view.get(c_e, p);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;
                        Scalar tr = 0.0;

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            el_energy += energy(this->params_, c[qp], c_grad_el[qp], el_strain.strain[qp], tr) * dx(qp);

                            if (this->params_.use_pressure) {
                                el_energy += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(
                                                 this->params_, c[qp]) *
                                             p[qp] * tr * dx(qp);
                            }
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x.comm().sum(val);

            assert(val == val);

            UTOPIA_TRACE_REGION_END("PhaseFieldForBrittleFractures::value(...)");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldForBrittleFractures::gradient(...)");

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
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *(this->local_pressure_field_));
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);
            auto press_val = press_fun.value(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            CoefStrain<USpace, Quadrature> strain(u_fun.coefficient(), q);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto p_view = press_val.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();
                auto ref_strain_u_view = ref_strain_u.view_device();

                auto g_view = this->space_.assembly_view_device(g);

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> stress_positive, stress_negative, stress;
                        StaticVector<Scalar, U_NDofs> u_el_vec;
                        StaticVector<Scalar, C_NDofs> c_el_vec;
                        Scalar tr_strain_u, gc, elast;

                        u_el_vec.set(0.0);
                        c_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);
                        auto u_grad_shape_el = v_grad_shape_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(c_e, c);

                        StaticVector<Scalar, NQuadPoints> p;
                        p_view.get(c_e, p);

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        for (int qp = 0; qp < NQuadPoints; ++qp) {
                            compute_stress(this->params_, c[qp], el_strain.strain[qp], stress, tr_strain_u, gc, elast);

                            for (SizeType j = 0; j < U_NDofs; ++j) {
                                auto &&strain_test = u_strain_shape_el(j, qp);
                                u_el_vec(j) += inner(stress, strain_test) * dx(qp);

                                if (this->params_.use_pressure) {
                                    u_el_vec(j) += gc * p[qp] * sum(diag(strain_test)) * dx(qp);
                                }
                            }

                            for (int j = 0; j < C_NDofs; ++j) {
                                const Scalar shape_test = c_shape_fun_el(j, qp);
                                const Scalar frac = PhaseFieldFracBase<FunctionSpace, Dim>::grad_fracture_energy_wrt_c(
                                    this->params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                                c_el_vec(j) += (elast * shape_test + frac) * dx(qp);

                                if (this->params_.use_pressure) {
                                    const Scalar der_c_pres =
                                        PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(
                                            this->params_, c[qp]) *
                                        p[qp] * tr_strain_u * shape_test;
                                    c_el_vec(j) += der_c_pres * dx(qp);
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

            this->space_.apply_zero_constraints(g);

            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(g, x_const);

                // // just a test...
                // auto* p_this =
                // const_cast<IsotropicPhaseFieldForBrittleFractures<FunctionSpace>
                // *>(this); add_irr_values_markers(p_this->_x_eq_values,
                // p_this->_eq_constrains_flg);
            }

            UTOPIA_TRACE_REGION_END("PhaseFieldForBrittleFractures::gradient(...)");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldForBrittleFractures::hessian(...)");

            if (empty(H)) {
                this->space_.create_matrix(H);
            } else {
                H *= 0.0;
            }

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            // udpate local pressure field
            this->space_.global_to_local(this->pressure_field_, *(this->local_pressure_field_));
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_pressure_field_);

            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto p_val = press_fun.value(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            CoefStrain<USpace, Quadrature> strain(u_fun.coefficient(), q);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = this->space_.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();
                auto p_view = p_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = this->space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> stress_positive;
                        StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;

                        MixedElem e;
                        space_view.elem(i, e);
                        el_mat.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                        auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> p;
                        c_view.get(c_e, c);
                        p_view.get(c_e, p);

                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        Scalar tr_strain_u, eep;

                        for (int qp = 0; qp < NQuadPoints; ++qp) {
                            compute_positive_quantities(
                                this->params_, el_strain.strain[qp], stress_positive, tr_strain_u, eep);

                            for (int l = 0; l < C_NDofs; ++l) {
                                const Scalar c_shape_l = c_shape_fun_el(l, qp);
                                auto &&c_grad_l = c_grad_shape_el(l, qp);
                                for (int j = l; j < C_NDofs; ++j) {
                                    Scalar val =
                                        PhaseFieldFracBase<FunctionSpace, Dim>::bilinear_cc(this->params_,
                                                                                            c[qp],
                                                                                            eep,
                                                                                            c_shape_fun_el(j, qp),
                                                                                            c_shape_l,
                                                                                            c_grad_shape_el(j, qp),
                                                                                            c_grad_l) *
                                        dx(qp);

                                    if (this->params_.use_pressure) {
                                        val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv2(
                                                   this->params_, c[qp]) *
                                               p[qp] * tr_strain_u * c_shape_fun_el(j, qp) * c_shape_l * dx(qp);
                                    }
                                    val = (l == j) ? (0.5 * val) : val;

                                    el_mat(l, j) += val;
                                    el_mat(j, l) += val;
                                }
                            }

                            for (SizeType l = 0; l < U_NDofs; ++l) {
                                auto &&u_strain_shape_l = u_strain_shape_el(l, qp);

                                for (SizeType j = l; j < U_NDofs; ++j) {
                                    Scalar val = bilinear_uu(this->params_,
                                                             c[qp],
                                                             el_strain.strain[qp],
                                                             u_strain_shape_el(j, qp),
                                                             u_strain_shape_l) *
                                                 dx(qp);
                                    val = (l == j) ? (0.5 * val) : val;
                                    el_mat(C_NDofs + l, C_NDofs + j) += val;
                                    el_mat(C_NDofs + j, C_NDofs + l) += val;
                                }
                            }

                            //////////////////////////////////////////////////////////////////////////////////////////////////////
                            if (this->params_.turn_off_uc_coupling == false ||
                                this->params_.turn_off_cu_coupling == false) {
                                for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                                    const Scalar c_shape_i = c_shape_fun_el(c_i, qp);
                                    for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                                        auto &&strain_shape = u_strain_shape_el(u_i, qp);

                                        Scalar val =
                                            PhaseFieldFracBase<FunctionSpace, Dim>::bilinear_uc(
                                                this->params_, c[qp], stress_positive, strain_shape, c_shape_i) *
                                            dx(qp);

                                        if (this->params_.use_pressure) {
                                            const Scalar tr_strain_shape = sum(diag(strain_shape));
                                            val += PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(
                                                       this->params_, c[qp]) *
                                                   p[qp] * tr_strain_shape * c_shape_i * dx(qp);
                                        }

                                        // not symetric, but more numerically stable... PF is weird thing to work
                                        // with....
                                        if (this->params_.turn_off_cu_coupling == false) {
                                            el_mat(c_i, C_NDofs + u_i) += val;
                                        }

                                        if (this->params_.turn_off_uc_coupling == false) {
                                            el_mat(C_NDofs + u_i, c_i) += val;
                                        }
                                    }
                                }
                            }
                        }

                        space_view.add_matrix(e, el_mat, H_view);
                    });
            }

            // check before boundary conditions
            // // if (this->check_derivatives_) {
            // this->diff_ctrl_.check_hessian(*this, x_const, H);
            // // }

            this->space_.apply_constraints(H);

            if (this->params_.use_crack_set_irreversibiblity) {
                this->apply_zero_constraints_irreversibiblity(H, x_const);
            }

            UTOPIA_TRACE_REGION_END("PhaseFieldForBrittleFractures::hessian(...)");
            return true;
        }

        //////////////////////////////////////////
        // for testing derivatives
        // void apply_zero_constraints_irreversibiblity(Matrix &H, const Vector &x) const override {
        //     std::vector<SizeType> indices;
        //     {
        //         Read<Vector> r(this->x_old_);
        //         Read<Vector> r2(x);

        //         Range range_w = range(this->x_old_);
        //         for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
        //             if (i % (Dim + 1) == 0) {
        //                 indices.push_back(i);
        //             }
        //         }
        //     }

        //     set_zero_rows(H, indices, 1.);
        // }

        // for testing derivatives
        // void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const override {
        //     {
        //         auto d_x_old = const_device_view(this->x_old_);
        //         auto d_x = const_device_view(x);

        //         auto g_view = view_device(g);
        //         parallel_for(
        //             range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
        //                 if (i % (Dim + 1) == 0) {
        //                     g_view.set(i, 0.0);
        //                 }
        //             });
        //     }
        // }

        template <class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(const PFFracParameters &params,
                                                    // c
                                                    const Scalar &phase_field_value,
                                                    const Grad &phase_field_grad,
                                                    // u
                                                    const Strain &strain,
                                                    Scalar &tr) {
            return PhaseFieldFracBase<FunctionSpace, Dim>::fracture_energy(
                       params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, strain, tr);
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(const PFFracParameters &params,
                                                            const Scalar &phase_field_value,
                                                            const Strain &strain,
                                                            Scalar &tr) {
            tr = trace(strain);
            Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Strain strain_dev2 = transpose(strain_dev) * strain_dev;
            Scalar tr2dv = trace(strain_dev2);

            // energy positive
            const Scalar energy_positive = (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * tr2dv);
            const Scalar energy_negative = (0.5 * params.kappa * tr_negative * tr_negative);

            const Scalar energy =
                (PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(params, phase_field_value) *
                 energy_positive) +
                energy_negative;

            return energy;
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const PFFracParameters &params,
                                                          const Scalar &phase_field_value,
                                                          const Strain &strain,
                                                          Stress &stress,
                                                          Scalar &tr,
                                                          Scalar &gc,
                                                          Scalar &elast_energy) {
            tr = trace(strain);
            const Strain strain_dev = strain - (((1. / Dim) * tr) * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            stress = ((params.kappa * tr_positive) * device::identity<Scalar>());
            stress += ((2.0 * params.mu) * strain_dev);

            gc = PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(params, phase_field_value);
            stress = (gc * stress) + ((params.kappa * tr_negative) * device::identity<Scalar>());

            // energy positive
            const Scalar energy_positive =
                (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));
            elast_energy =
                PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation_deriv(params, phase_field_value) *
                energy_positive;
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Grad &strain,
                                                         const Grad &strain_trial,
                                                         const Grad &strain_test) {
            const Scalar strain0tr = trace(strain);

            Tensor4th<Scalar, Dim, Dim, Dim, Dim> Jacobian_neg, Jacobian_mult;

            if (strain0tr < 0) {
                Jacobian_neg = params.kappa * params.I4sym;
            }

            Scalar gc = PhaseFieldFracBase<FunctionSpace, Dim>::quadratic_degradation(params, phase_field_value);

            // would be nicer, if this works without 4th order tensor...
            Jacobian_mult = params.elast_tensor - Jacobian_neg;
            Jacobian_mult = (gc * Jacobian_mult) + Jacobian_neg;

            Scalar val = inner(strain_trial, contraction(Jacobian_mult, strain_test));

            return val;
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_positive_quantities(const PFFracParameters &params,
                                                                       const Strain &strain,
                                                                       Stress &stress_positive,
                                                                       Scalar &tr,
                                                                       Scalar &energy_positive) {
            tr = trace(strain);
            const Strain strain_dev = strain - (((1. / Dim) * tr) * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            stress_positive = ((params.kappa * tr_positive) * device::identity<Scalar>());
            stress_positive += ((2.0 * params.mu) * strain_dev);

            energy_positive =
                (0.5 * params.kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));
        }
    };

}  // namespace utopia
#endif