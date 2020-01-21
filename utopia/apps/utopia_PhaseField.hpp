#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    template<class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PhaseFieldForBrittleFractures final : public Function<
                            typename FunctionSpace::Matrix,
                            typename FunctionSpace::Vector>

    {
    public:
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector   = typename FunctionSpace::Vector;
        using Matrix   = typename FunctionSpace::Matrix;
        using Device   = typename FunctionSpace::Device;

        using USpace   = typename FunctionSpace::template Subspace<Dim>;
        using CSpace   = typename FunctionSpace::template Subspace<1>;

        using UElem    = typename USpace::ViewDevice::Elem;
        using CElem    = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;

        //FIXME
        using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        class Parameters {
        public:

            Parameters()
            : a(1.0), b(1.0), d(1.0), f(1.0), length_scale(1.0), fracture_toughness(1.0), mu(1.0), lambda(1.0)
            {}

            Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda;

        };

        PhaseFieldForBrittleFractures(FunctionSpace &space)
        : space_(space)
        {}

        // void assemble(
        //     Vector &x,
        //     Matrix &H,
        //     Vector &g,
        //     Scalar &val
        // )
        // {
        //     USpace U;
        //     space_.subspace(0, U);
        //     CSpace C = space_.subspace(Dim);





        // }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const
        {
            space_.create_matrix(H);
            return true;
        }

        inline bool update(const Vector &x) override {
            // x_coeff_.update(x);
            return true;
        }

        bool value(
            const Vector &x_const,
            Scalar &val
        ) const override
        {
            USpace U;
            space_.subspace(0, U);
            CSpace C = space_.subspace(Dim);

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val  = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            PrincipalStrains<USpace, Quadrature> strain(U, q);
            strain.update(x);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();

                auto c_view      = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view      = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
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

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar sum_eigs = sum(el_strain.values[qp]);
                            strain_view.split(el_strain, qp, strain_n, strain_p);

                            el_energy += energy(params_, c[qp], c_grad_el[qp], sum_eigs, strain_n, strain_p) * dx(qp);

                        }

                        return el_energy;
                    },
                    val
                );
            }

            val = x.comm().sum(val);

            disp(val);
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override
        {
            if(empty(g)) {
                space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            USpace U;
            space_.subspace(0, U);
            CSpace C = space_.subspace(Dim);

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val  = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape      = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(U, q);
            strain.update(x);

            {
                auto U_view      = U.view_device();
                auto C_view      = C.view_device();

                auto c_view      = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view      = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = space_.assembly_view_device(g);

                Device::parallel_for(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
                        StaticMatrix<Scalar, Dim, Dim> stress, strain_p;
                        StaticVector<Scalar, U_NDofs> u_el_vec;
                        StaticVector<Scalar, C_NDofs>  c_el_vec;

                        u_el_vec.set(0.0);
                        c_el_vec.set(0.0);

                        ////////////////////////////////////////////

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto u_grad_shape_el = v_grad_shape_view.make(u_e);

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(c_e, c);

                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx        = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el  = c_shape_view.make(c_e);

                        ////////////////////////////////////////////

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar sum_eigs = sum(el_strain.values[qp]);

                            strain_view.split_positive(el_strain, qp, strain_p);

                            split_stress(params_, c[qp], el_strain.values[qp], el_strain.vectors[qp], stress);

                            for(SizeType j = 0; j < u_grad_shape_el.n_functions(); ++j) {
                                u_el_vec(j) += inner(stress, u_grad_shape_el(j, qp)) * dx(qp);
                            }

                            const Scalar elast =
                                grad_elastic_energy_wrt_c(
                                            params_,
                                            c[qp],
                                            c_grad_el[qp],
                                            sum_eigs,
                                            strain_p
                                );


                            for(SizeType j = 0; j < c_grad_shape_el.n_functions(); ++j) {
                                const Scalar tf   = c_shape_fun_el(j, qp);
                                const Scalar frac =
                                    grad_fracture_energy_wrt_c(
                                            params_,
                                            c[qp],
                                            c_grad_el[qp],
                                            tf,
                                            c_grad_shape_el(j, qp)
                                    );

                                c_el_vec(j) += (elast * tf + frac) * dx(qp);
                            }

                        }

                        U_view.add_vector(u_e, u_el_vec, g_view);
                        C_view.add_vector(c_e, c_el_vec, g_view);
                    }
                );
            }

            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override
        {
            if(empty(H)) {
                space_.create_matrix(H);
            } else {
                H *= 0.0;
            }


            USpace U;
            space_.subspace(0, U);
            CSpace C = space_.subspace(Dim);

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val  = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape      = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(U, q);
            strain.update(x);

            {
                auto U_view      = U.view_device();
                auto C_view      = C.view_device();
                auto space_view  = space_.view_device();

                auto c_view      = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view      = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = space_.assembly_view_device(H);

                Device::parallel_for(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
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

                        ////////////////////////////////////////////

                        CElem c_e;
                        C_view.elem(i, c_e);
                        StaticVector<Scalar, NQuadPoints> c;
                        c_view.get(c_e, c);

                        auto dx        = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el  = c_shape_view.make(c_e);

                        ////////////////////////////////////////////

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar sum_eigs = sum(el_strain.values[qp]);
                            strain_view.split(el_strain, qp, strain_n, strain_p);

                            const Scalar eep = elastic_energy_positve(params_, sum_eigs, strain_p);

                            for(SizeType l = 0; l < c_grad_shape_el.n_functions(); ++l) {
                                for(SizeType j = 0; j < c_grad_shape_el.n_functions(); ++j) {
                                    el_mat(l, j) += (
                                        diffusion_c(params_, c_grad_shape_el(j, qp), c_grad_shape_el(l, qp)) +
                                        reaction_c(params_,  c_shape_fun_el(j, qp),  c_shape_fun_el(l, qp))  +
                                        elastic_deriv_cc(params_, c[qp], eep, c_shape_fun_el(l, qp))
                                    ) * dx(qp);
                                }
                            }

                            for(SizeType l = 0; l < u_grad_shape_el.n_functions(); ++l) {
                                for(SizeType j = 0; j < u_grad_shape_el.n_functions(); ++j) {
                                    el_mat(C_NDofs + l, C_NDofs + j) += bilinear_uu(
                                        params_,
                                        c[qp],
                                        sum_eigs,
                                        strain_n,
                                        strain_p,
                                        u_grad_shape_el(j, qp),
                                        u_grad_shape_el(l, qp)
                                        ) * dx(qp);
                                }
                            }

                            for(SizeType u_i = 0; u_i < u_grad_shape_el.n_functions(); ++u_i) {
                                for(SizeType c_i = 0; c_i < c_grad_shape_el.n_functions(); ++c_i) {

                                    const Scalar val = bilinear_cu(
                                                params_,
                                                c[qp],
                                                sum_eigs,
                                                strain_p,
                                                strain_p - strain_n,
                                                c_shape_fun_el(c_i, qp),
                                                u_grad_shape_el(u_i, qp)
                                            ) * dx(qp);


                                    el_mat(c_i, C_NDofs + u_i) += val;
                                    el_mat(C_NDofs + u_i, c_i) += val;
                                }
                            }
                        }

                        space_view.add_matrix(e, el_mat, H_view);
                    }
                );
            }


            return true;
        }

        //////////////////////////////////////////

        template<class Strain, class FullStrain, class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cu(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &trace,
            const Strain &strain_positive,
            const FullStrain &full_strain,
            const Scalar &c_trial_fun,
            const Grad &u_grad_test
            )
        {
            auto C_test  = 0.5 * (u_grad_test  + transpose(u_grad_test));

            const Scalar trace_positive = split_p(trace);
            const Scalar trial =
                strain_energy(params, trace_positive, strain_positive) *
                quadratic_degradation_deriv(
                           params,
                            phase_field_value) * c_trial_fun;

            return trial * inner(full_strain, C_test);

        }

        template<class Strain, class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &trace,
            const Strain &strain_negative,
            const Strain &strain_positive,
            const Grad &g_trial,
            const Grad &g_test
            )
        {
            auto C_trial = 0.5 * (g_trial + transpose(g_trial));
            auto C_test  = 0.5 * (g_test  + transpose(g_test));

            const Scalar contr_C = inner(C_trial, C_test);

            const Scalar trace_negative = split_n(trace);
            const Scalar trace_positive = split_p(trace);

            const Scalar positive_part = strain_energy(params, trace_positive, strain_positive) * contr_C;
            const Scalar negative_part = strain_energy(params, trace_negative, strain_negative) * contr_C;
            return quadratic_degradation(params, phase_field_value) * positive_part + negative_part;
        }

        template<class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar diffusion_c(
            const Parameters &params,
            const Grad &g_trial,
            const Grad &g_test
            )
        {
            return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
        }

        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(
            const Parameters &params,
            const Scalar &trial,
            const Scalar &test
            )
        {
            return (params.fracture_toughness / params.length_scale) * trial * test;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &elastic_energy_positive,
            const Scalar &test
            )
        {
            const Scalar dcc = quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * elastic_energy_positive * test;
        }


        template<class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar grad_elastic_energy_wrt_c(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad,
            const Scalar &trace,
            const Strain &strain_positive
            )
        {
            const Scalar trace_positive = split_p(trace);
            return quadratic_degradation_deriv(
                params, phase_field_value) * strain_energy(params, trace_positive, strain_positive);
        }

        template<class Grad, class GradTest>
        UTOPIA_INLINE_FUNCTION static Scalar grad_fracture_energy_wrt_c(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad,
            const Scalar &test_function,
            const GradTest &grad_test_function
            )
        {
            return params.fracture_toughness * (
                1./params.length_scale * phase_field_value * test_function +
                params.length_scale * inner(phase_field_grad, grad_test_function)
            );
        }

        template<class EigenValues, class EigenMatrix, class Stress>
        UTOPIA_INLINE_FUNCTION static void split_stress(
            const Parameters &params,
            const Scalar &phase_field_value,
            const EigenValues &values,
            const EigenMatrix &mat,
            Stress &stress
            )
        {
            Scalar tr = sum(values);
            const Scalar tr_p = split_p(tr);
            const Scalar tr_n = split_n(tr);

            StaticVector<Scalar, Dim> v;

            stress.set(0.0);

            for(int d = 0; d < Dim; ++d) {
                const Scalar eig_p = split_p(values[d]);
                const Scalar eig_n = split_n(values[d]);

                const Scalar val_p = quadratic_degradation(params, phase_field_value) * (params.lambda * tr_p + 2.0 * params.mu * eig_p);
                const Scalar val_n = params.lambda * tr_n + 2.0 * params.mu * eig_n;

                const Scalar val = val_p + val_n;

                mat.col(d, v);
                stress += val * outer(v, v);
            }
        }

        template<class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(
            const Parameters &params,
            //c
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad,
            // u
            const Scalar &trace,
            const Strain &strain_n,
            const Strain &strain_p)
        {
            return
            fracture_energy(params, phase_field_value, phase_field_grad) +
            elastic_energy(params,  phase_field_value, trace, strain_n, strain_p);
        }

        UTOPIA_INLINE_FUNCTION static Scalar split_p(const Scalar &x) {
            return (device::abs(x) + x)/2;
        };

        UTOPIA_INLINE_FUNCTION static Scalar split_n(const Scalar &x) {
            return (device::abs(x) - x)/2;
        };

        template<class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar fracture_energy(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad
            )
        {
            return params.fracture_toughness * (
                1./(2.0 * params.length_scale) * phase_field_value * phase_field_value +
                params.length_scale/2.0 * inner(phase_field_grad, phase_field_grad)
            );
        }

        template<class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar strain_energy(
            const Parameters &params,
            const Scalar trace,
            const Strain &strain)
        {
            Scalar tr = trace;
            return 0.5 * params.lambda * tr*tr + params.mu * inner(strain, strain);
        }

        template<class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &trace,
            const Strain &strain_negative,
            const Strain &strain_positive
            )
        {
            const Scalar trace_negative = split_n(trace);
            const Scalar trace_positive = split_p(trace);
            return
                quadratic_degradation(params, phase_field_value) * strain_energy(params, trace_positive, strain_positive) +
                strain_energy(params, trace_negative, strain_negative);
        }

        template<class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy_positve(
            const Parameters &params,
            const Scalar &trace,
            const Strain &strain_positive
            )
        {
            const Scalar trace_positive = split_p(trace);
            return strain_energy(params, trace_positive, strain_positive);
        }

        UTOPIA_INLINE_FUNCTION static Scalar degradation(
            const Parameters &params,
            const Scalar &c
            )
        {
            Scalar imc = 1.0 - c;
            Scalar res = params.f + imc * params.d;
            imc *= imc;
            res += params.b * imc;
            imc *= imc; //FIXME

            res += params.a * imc;
            return res;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation(
           const Parameters &,
            const Scalar &c
            )
        {
            Scalar imc = 1.0 - c;
            return imc*imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv(
           const Parameters &,
            const Scalar &c
            )
        {
            Scalar imc = 1.0 - c;
            return -2.0 * imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv2(
           const Parameters &,
            const Scalar &
            )
        {
            return 2.0;
        }

    private:
        FunctionSpace space_;
        Parameters params_;
    };

    // template<class FunctionSpace>
    // static void compute_strain_energy_splitting(
    //     const FunctionSpace &space,
    //     const typename FunctionSpace::Vector &displacement)
    // {
    //     static const int Dim = FunctionSpace::Dim;
    //     static const int NVars = Dim;

    //     using Elem           = typename FunctionSpace::Elem;
    //     using ElemView       = typename FunctionSpace::ViewDevice::Elem;
    //     using Mesh           = typename FunctionSpace::Mesh;
    //     using SizeType       = typename Mesh::SizeType;
    //     using Scalar         = typename Mesh::Scalar;
    //     using Quadrature     = utopia::Quadrature<Elem, 2>;
    //     using Dev            = typename FunctionSpace::Device;
    //     using VectorD        = utopia::StaticVector<Scalar, Dim>;
    //     using MatrixDxD      = utopia::StaticMatrix<Scalar, Dim, Dim>;

    //     Quadrature q;
    //     PrincipalStrains<FunctionSpace, Quadrature> strain(space, q);
    //     strain.update(displacement);

    //     Differential<FunctionSpace, Quadrature> differential(space, q);

    //     Scalar mu = 1.0, lambda = 1.0;



    //     //end of host code

    //     StaticVector<Scalar, 2> energy;
    //     energy.set(0.0);
    //     {
    //         //beginning of device code
    //         auto strain_view = strain.view_device();
    //         auto space_view  = space.view_device();
    //         auto differential_view = differential.view_device();

    //         Dev::parallel_reduce(
    //             space.local_element_range(),
    //             UTOPIA_LAMBDA(const SizeType &i)
    //         {
    //             VectorD v;
    //             StaticVector<Scalar, 2> e_energy;
    //             MatrixDxD strain_n;
    //             MatrixDxD strain_p;

    //             ElemView e;
    //             space_view.elem(i, e);

    //             auto el_strain = strain_view.make(e);
    //             auto dx        = differential_view.make(e);

    //             const SizeType n_qp = el_strain.values.size();

    //             e_energy.set(0.0);

    //             for(SizeType qp = 0; qp < n_qp; ++qp) {
    //                 //reset strain tensor
    //                 strain_n.set(0.0);
    //                 strain_p.set(0.0);

    //                 //compute splitted strain
    //                 for(int d = 0; d < Dim; ++d) {
    //                     auto e_val = el_strain.values[qp][d];
    //                     el_strain.vectors[qp].col(d, v);

    //                     auto outer_v = outer(v, v);

    //                     auto eig_p = split_p(e_val);
    //                     auto eig_n = split_m(e_val);

    //                     strain_n += eig_n * outer_v;
    //                     strain_p += eig_p * outer_v;
    //                 }

    //                 const Scalar trace_e_n = trace(strain_n);
    //                 const Scalar trace_e_p = trace(strain_p);

    //                 e_energy[0] += ((0.5 * lambda * trace_e_n * trace_e_n) +
    //                            (mu * inner(strain_n, strain_n))) * dx(qp);

    //                 e_energy[1] += ((0.5 * lambda * trace_e_p * trace_e_p) +
    //                            (mu * inner(strain_p, strain_p))) * dx(qp);
    //             }

    //             return e_energy;
    //         }, energy);
    //     }

    //     space.comm().sum(2, &energy[0]);

    //     if(space.comm().rank() == 0) {
    //         disp(energy);
    //     }
    // }
}
#endif