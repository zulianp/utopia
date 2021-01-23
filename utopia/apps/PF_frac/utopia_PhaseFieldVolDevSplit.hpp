#ifndef UTOPIA_PHASE_FIELD_VOL_DEV_SPLIT_HPP
#define UTOPIA_PHASE_FIELD_VOL_DEV_SPLIT_HPP

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
        // using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;
        using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 0>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        PhaseFieldVolDevSplit(FunctionSpace &space) : PhaseFieldFracBase<FunctionSpace, Dim>(space) {}

        PhaseFieldVolDevSplit(FunctionSpace &space, const PFFracParameters &params)
            : PhaseFieldFracBase<FunctionSpace, Dim>(space, params) {}

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("PhaseFieldVolDevSplit::value(...)");

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;
            PrincipalStrains<USpace, Quadrature> strain(u_fun.coefficient(), q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

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
                        c_view.get(c_e, c);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;

                        // std::cout << "NQuadPoints: " << NQuadPoints << "   \n";
                        // exit(0);

                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            // el_energy += energy(this->params_, c[qp], c_grad_el[qp], el_strain.strain[qp]) * dx(qp);
                            el_energy += elastic_energy(this->params_, c[qp], el_strain.strain[qp]) * dx(qp);
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

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(u_fun.coefficient(), q);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

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
                        StaticMatrix<Scalar, Dim, Dim> stress_positive, stress_negative;
                        StaticVector<Scalar, U_NDofs> u_el_vec;
                        StaticVector<Scalar, C_NDofs> c_el_vec;

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

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////

                        // std::cout << "NQuadPoints:  " << NQuadPoints << " \n";
                        for (int qp = 0; qp < NQuadPoints; ++qp) {
                            compute_stress(
                                this->params_, c[qp], el_strain.strain[qp], stress_positive, stress_negative);

                            for (SizeType j = 0; j < U_NDofs; ++j) {
                                auto &&strain_test = u_strain_shape_el(j, qp);

                                u_el_vec(j) += inner(stress_positive, strain_test) * dx(qp);
                                u_el_vec(j) += inner(stress_negative, strain_test) * dx(qp);
                            }

                            // exit(0);

                            // const Scalar elast = grad_elastic_energy_wrt_c(this->params_, c[qp],
                            // el_strain.strain[qp]);

                            // for (int j = 0; j < C_NDofs; ++j) {
                            //     const Scalar shape_test = c_shape_fun_el(j, qp);
                            //     const Scalar frac = grad_fracture_energy_wrt_c(
                            //         this->params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                            //     c_el_vec(j) += (elast * shape_test + frac) * dx(qp);
                            // }
                        }

                        U_view.add_vector(u_e, u_el_vec, g_view);
                        // C_view.add_vector(c_e, c_el_vec, g_view);
                    });
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_grad(*this, x_const, g);
            }

            this->space_.apply_zero_constraints(g);

            if (this->params_.use_crack_set_irreversibiblity) {
                apply_zero_constraints_irreversibiblity(g, x_const);

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

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(u_fun.coefficient(), q);
            PrincipalShapeStress<USpace, Quadrature> p_stress(U, q, this->params_.mu, this->params_.lambda);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = this->space_.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                // FIXME
                auto p_stress_view = p_stress.view_device();

                auto H_view = this->space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    this->space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        StaticMatrix<Scalar, Dim, Dim> strain_n, strain_p;
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
                        c_view.get(c_e, c);

                        auto dx = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        const int n_u_fun = u_grad_shape_el.n_functions();
                        const int n_c_fun = c_grad_shape_el.n_functions();

                        for (int qp = 0; qp < NQuadPoints; ++qp) {
                            // const Scalar eep = elastic_energy_positive(this->params_, el_strain.strain[qp]);

                            // for (int l = 0; l < n_c_fun; ++l) {
                            //     for (int j = 0; j < n_c_fun; ++j) {
                            //         el_mat(l, j) += bilinear_cc(this->params_,
                            //                                     c[qp],
                            //                                     eep,
                            //                                     c_shape_fun_el(j, qp),
                            //                                     c_shape_fun_el(l, qp),
                            //                                     c_grad_shape_el(j, qp),
                            //                                     c_grad_shape_el(l, qp)) *
                            //                         dx(qp);
                            //     }
                            // }

                            //
                            for (int l = 0; l < U_NDofs; ++l) {
                                auto &&u_strain_shape_l = u_strain_shape_el(l, qp);
                                for (int j = 0; j < U_NDofs; ++j) {
                                    el_mat(C_NDofs + l, C_NDofs + j) +=
                                        bilinear_uu(this->params_, c[qp], u_strain_shape_el(j, qp), u_strain_shape_l) *
                                        dx(qp);
                                }
                            }

                            //////////////////////////////////////////////////////////////////////////////////////////////////////
                            // if (this->params_.turn_off_uc_coupling == false ||
                            //     this->params_.turn_off_cu_coupling == false) {
                            //     StaticMatrix<Scalar, Dim, Dim> stress_positive;
                            //     compute_stress_positive(this->params_, el_strain.strain[qp], stress_positive);

                            //     for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                            //         const Scalar c_shape_i = c_shape_fun_el(c_i, qp);
                            //         for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                            //             auto &&strain_shape = u_strain_shape_el(u_i, qp);

                            //             Scalar val =
                            //                 bilinear_uc(
                            //                     this->params_, c[qp], stress_positive, strain_shape, c_shape_i) *
                            //                 dx(qp);

                            //             // not symetric, but more numerically stable
                            //             if (this->params_.turn_off_cu_coupling == false) {
                            //                 el_mat(c_i, C_NDofs + u_i) += val;
                            //             }

                            //             if (this->params_.turn_off_uc_coupling == false) {
                            //                 el_mat(C_NDofs + u_i, c_i) += val;
                            //             }
                            //         }
                            //     }
                            // }
                        }

                        space_view.add_matrix(e, el_mat, H_view);
                    });
            }

            // check before boundary conditions
            if (this->check_derivatives_) {
                this->diff_ctrl_.check_hessian(*this, x_const, H);
            }

            this->space_.apply_constraints(H);

            if (this->params_.use_crack_set_irreversibiblity) {
                apply_zero_constraints_irreversibiblity(H, x_const);
            }

            UTOPIA_TRACE_REGION_END("PhaseFieldForBrittleFractures::hessian(...)");
            return true;
        }

        //////////////////////////////////////////

        void apply_zero_constraints_irreversibiblity(Matrix &H, const Vector &x) const override {
            std::vector<SizeType> indices;
            {
                Read<Vector> r(this->x_old_);
                Read<Vector> r2(x);

                Range range_w = range(this->x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 0) {
                        indices.push_back(i);
                    }
                }
            }

            set_zero_rows(H, indices, 1.);
        }

        void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const override {
            {
                auto d_x_old = const_device_view(this->x_old_);
                auto d_x = const_device_view(x);

                auto g_view = view_device(g);
                parallel_for(
                    range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            g_view.set(i, 0.0);
                        }
                    });
            }
        }

        template <class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Scalar &elastic_energy_p,
                                                         const Scalar &shape_trial,
                                                         const Scalar &shape_test,
                                                         const GradShape &grad_trial,
                                                         const GradShape &grad_test) {
            return diffusion_c(params, grad_trial, grad_test) + reaction_c(params, shape_trial, shape_test) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy_p, shape_trial, shape_test);
        }

        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return quadratic_degradation_deriv(params, phase_field_value) * c_trial_fun * inner(stress_p, full_strain);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar diffusion_c(const PFFracParameters &params,
                                                         const Grad &g_trial,
                                                         const Grad &g_test) {
            return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
        }

        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(const PFFracParameters &params,
                                                        const Scalar &trial,
                                                        const Scalar &test) {
            return (params.fracture_toughness / params.length_scale) * trial * test;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const PFFracParameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &elastic_energy_positive,
                                                              const Scalar &trial,
                                                              const Scalar &test) {
            const Scalar dcc = quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * trial * elastic_energy_positive * test;
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar grad_elastic_energy_wrt_c(const PFFracParameters &params,
                                                                       const Scalar &phase_field_value,
                                                                       const Strain &strain) {
            const Scalar tr = trace(strain);
            const Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Scalar kappa = params.lambda + (2.0 * params.mu / Dim);

            // energy positive
            const Scalar energy_positive =
                (0.5 * kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));

            return quadratic_degradation_deriv(params, phase_field_value) * energy_positive;
        }

        template <class Grad, class GradTest>
        UTOPIA_INLINE_FUNCTION static Scalar grad_fracture_energy_wrt_c(const PFFracParameters &params,
                                                                        const Scalar &phase_field_value,
                                                                        const Grad &phase_field_grad,
                                                                        const Scalar &test_function,
                                                                        const GradTest &grad_test_function) {
            return params.fracture_toughness * (1. / params.length_scale * phase_field_value * test_function +
                                                params.length_scale * inner(phase_field_grad, grad_test_function));
        }

        template <class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(const PFFracParameters &params,
                                                    // c
                                                    const Scalar &phase_field_value,
                                                    const Grad &phase_field_grad,
                                                    // u
                                                    const Strain &strain) {
            return fracture_energy(params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, strain);
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(const PFFracParameters &params,
                                                            const Scalar &phase_field_value,
                                                            const Strain &strain) {
            const Scalar tr = trace(strain);
            Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Scalar kappa = params.lambda + (2.0 * params.mu / Dim);
            const Strain strain_dev2 = transpose(strain_dev) * strain_dev;
            Scalar tr2dv = trace(strain_dev2);

            // energy positive
            const Scalar energy_positive = (0.5 * kappa * tr_positive * tr_positive) + (params.mu * tr2dv);
            const Scalar energy_negative = (0.5 * kappa * tr_negative * tr_negative);

            const Scalar energy =
                (quadratic_degradation(params, phase_field_value) * energy_positive) + energy_negative;

            return energy;
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy_positive(const PFFracParameters &params,
                                                                     const Strain &strain) {
            const Scalar tr = trace(strain);
            const Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Scalar kappa = params.lambda + (2.0 * params.mu / Dim);

            // energy positive
            const Scalar energy_positive =
                (0.5 * kappa * tr_positive * tr_positive) + (params.mu * inner(strain_dev, strain_dev));

            return energy_positive;
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const PFFracParameters &params,
                                                          const Scalar &phase_field_value,
                                                          const Strain &strain,
                                                          Stress &stress_positive,
                                                          Stress &stress_negative) {
            const Scalar tr = trace(strain);
            const Strain strain_dev = strain - (((1. / Dim) * tr) * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Scalar kappa = params.lambda + ((2.0 * params.mu) / Dim);

            stress_positive = ((kappa * tr_positive) * device::identity<Scalar>());
            stress_positive += ((2.0 * params.mu) * strain_dev);
            // stress_positive = quadratic_degradation(params, phase_field_value) * stress_positive;

            stress_negative = ((kappa * tr_negative) * device::identity<Scalar>());

            // std::cout << "kappa: " << kappa << "  \n";
            // disp(stress_positive);
            // std::cout << " \n \n \n";
            // disp(stress_negative);
            // exit(0);

            // stress = (quadratic_degradation(params, phase_field_value) * stress_positive) + stress_negative;
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Grad &strain,
                                                         const Grad &strain_test) {
            // const StaticMatrix<Scalar, Dim, Dim> strain = 0.5 * (g_trial + transpose(g_trial));

            StaticMatrix<Scalar, Dim, Dim> stress_positive, stress_negative;
            compute_stress(params, phase_field_value, strain, stress_positive, stress_negative);

            // const Scalar gc = quadratic_degradation(params, phase_field_value);

            return (inner(stress_positive, strain_test)) + inner(stress_negative, strain_test);
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress_positive(const PFFracParameters &params,
                                                                   const Strain &strain,
                                                                   Stress &stress_positive) {
            const Scalar tr = trace(strain);
            const Strain strain_dev = strain - ((1. / Dim) * tr * device::identity<Scalar>());

            const Scalar tr_negative = device::min(tr, 0.0);
            const Scalar tr_positive = tr - tr_negative;

            const Scalar kappa = params.lambda + (2.0 * params.mu / Dim);

            stress_positive = (kappa * tr_positive * device::identity<Scalar>());
            stress_positive += (2.0 * params.mu * strain_dev);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar fracture_energy(const PFFracParameters &params,
                                                             const Scalar &phase_field_value,
                                                             const Grad &phase_field_grad) {
            return params.fracture_toughness *
                   (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
                    params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation(const PFFracParameters &, const Scalar &c) {
            Scalar imc = 1.0 - c;
            return imc * imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv(const PFFracParameters &, const Scalar &c) {
            Scalar imc = 1.0 - c;
            return -2.0 * imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv2(const PFFracParameters &, const Scalar &) {
            return 2.0;
        }
    };

}  // namespace utopia
#endif