#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_Views.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_TensorView4.hpp"
#include "utopia_DeviceTensorProduct.hpp"
#include "utopia_DeviceTensorContraction.hpp"
#include "utopia_StrainView.hpp"

namespace utopia {

    template<class FunctionSpace, int Dim = FunctionSpace::Dim>
    class IsotropicPhaseFieldForBrittleFractures final : public Function<
                            typename FunctionSpace::Matrix,
                            typename FunctionSpace::Vector>,
                            public Configurable

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

        class Parameters : public Configurable {
        public:

            void read(Input &in) override
            {
                in.get("a", a);
                in.get("b", b);
                in.get("d", d);
                in.get("f", f);
                in.get("length_scale", length_scale);
                in.get("fracture_toughness", fracture_toughness);
                in.get("mu", mu);
                in.get("lambda", lambda);
                in.get("regularization", regularization);
                in.get("pressure", pressure);

                in.get("use_penalty_irreversibility", use_penalty_irreversibility);
                in.get("penalty_param", penalty_param);

                in.get("use_crack_set_irreversibiblity", use_crack_set_irreversibiblity); 
                in.get("crack_set_tol", crack_set_tol); 
            }


            Parameters()
            :   a(1.0), b(1.0), d(1.0), f(1.0), length_scale(1.0), fracture_toughness(1.0), 
                mu(1.0), lambda(1.0), regularization(1e-10), pressure(0.0), penalty_param(0.0), 
                crack_set_tol(0.95), use_penalty_irreversibility(false), use_crack_set_irreversibiblity(false)
            {}

            Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda; 
            Scalar regularization, pressure, penalty_param, crack_set_tol;
            bool use_penalty_irreversibility, use_crack_set_irreversibiblity; 
        };

        void read(Input &in) override
        {
            params_.read(in);
            in.get("use_dense_hessian", use_dense_hessian_);
            in.get("check_derivatives", check_derivatives_);
            in.get("diff_controller", diff_ctrl_);
        }

        // this is a bit of hack 
        void set_pressure(const Scalar & pressure){
            params_.pressure = pressure; 
        }


        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space)
        : space_(space), use_dense_hessian_(false), check_derivatives_(false)
        {
            params_.length_scale = 2.0 * space.mesh().min_spacing();
            params_.fracture_toughness = 0.001;

            params_.mu = 80.0;
            params_.lambda = 120.0;

            // this computation follows eq. 50 from "On penalization in variational phase-field models of britlle fracture, Gerasimov, Lorenzis"
            if(params_.use_penalty_irreversibility)
            {
                Scalar tol = 1e-3; 
                Scalar tol2 = tol*tol; 
                params_.penalty_param = params_.fracture_toughness/params_.length_scale * (1.0/ tol2 - 1.0);             
            }
        }

        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space, const Parameters &params)
        : space_(space), params_(params), use_dense_hessian_(false), check_derivatives_(false)
        {}

        void use_dense_hessian(const bool val)
        {
            use_dense_hessian_ = val;
        }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const
        {
            if(use_dense_hessian_) {
                H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
            } else {
                space_.create_matrix(H);
            }
            return true;
        }

        inline bool update(const Vector &x) override {
            // x_coeff_.update(x);
            return true;
        }

        bool value(const Vector &x_const, Scalar &val) const override
        {
            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            auto &x_old = const_cast<Vector &>(x_old_);
            FEFunction<CSpace> c_old_fun(C, x_old);            

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_old  = c_old_fun.value(q);             

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
                auto c_old_view  = c_old.view_device(); 

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
                        StaticVector<Scalar, NQuadPoints> c_old;
                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old); 

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {

                            Scalar tr = trace(el_strain.strain[qp]); 
                            if(params_.pressure > 0){
                                el_energy += quadratic_degradation(params_,  c[qp]) *  params_.pressure *  tr * dx(qp); 
                            }


                            el_energy += energy(params_, c[qp], c_grad_el[qp], tr, el_strain.strain[qp]) * dx(qp);

                            if(params_.use_penalty_irreversibility){
                                auto c_cold             = c[qp] - c_old[qp]; 
                                auto c_cold_bracket     = c_cold < 0.0 ? c_cold: 0.0; 
                                el_energy               += params_.penalty_param/2.0 * c_cold_bracket * c_cold_bracket * dx(qp);
                            }                            
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val
                );
            }

            val = x.comm().sum(val);

            assert(val == val);
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
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            auto &x_old = const_cast<Vector &>(x_old_);
            FEFunction<CSpace> c_old_fun(C, x_old);            

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_old  = c_old_fun.value(q); 

            auto c_grad = c_fun.gradient(q);
            auto u_val  = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape      = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(U, q);
            strain.update(x);

            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view      = U.view_device();
                auto C_view      = C.view_device();

                auto c_view      = c_val.view_device();
                auto c_old_view  = c_old.view_device();   

                auto c_grad_view = c_grad.view_device();
                auto u_view      = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = space_.assembly_view_device(g);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
                        StaticMatrix<Scalar, Dim, Dim> stress, strain_p;
                        StaticVector<Scalar, U_NDofs> u_el_vec;
                        StaticVector<Scalar, C_NDofs>  c_el_vec;

                        // StaticMatrix<Scalar, Dim, Dim> identity; 
                        // identity.identity();                         

                        u_el_vec.set(0.0);
                        c_el_vec.set(0.0);

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

                        StaticVector<Scalar, NQuadPoints> c_old;
                        c_old_view.get(c_e, c_old);


                        auto c_grad_el = c_grad_view.make(c_e);
                        auto dx        = differential_view.make(c_e);
                        auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el  = c_shape_view.make(c_e);

                        ////////////////////////////////////////////

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {

                            const Scalar tr_strain_u = trace(el_strain.strain[qp]); 

                            compute_stress(params_, tr_strain_u, el_strain.strain[qp],  stress); 
                            stress = (quadratic_degradation(params_, c[qp]) * (1.0 - params_.regularization) + params_.regularization) * stress; 

                            for(SizeType j = 0; j < U_NDofs; ++j) {
                                auto &&strain_test = u_strain_shape_el(j, qp);
                                u_el_vec(j) += inner(stress, strain_test) * dx(qp);

                                if(params_.pressure > 0){
                                    u_el_vec(j) +=  quadratic_degradation(params_, c[qp]) *  params_.pressure  * sum(diag(strain_test)) * dx(qp);
                                }
                            }


                            const Scalar elast =
                                grad_elastic_energy_wrt_c(
                                            params_,
                                            c[qp],
                                            tr_strain_u,
                                            el_strain.strain[qp]
                                );


                            for(SizeType j = 0; j < C_NDofs; ++j) {
                                const Scalar shape_test   = c_shape_fun_el(j, qp);
                                const Scalar frac =
                                    grad_fracture_energy_wrt_c(
                                            params_,
                                            c[qp],
                                            c_grad_el[qp],
                                            shape_test,
                                            c_grad_shape_el(j, qp)
                                    );

                                
                                if(params_.pressure > 0){
                                    const Scalar der_c_pres  = quadratic_degradation_deriv(params_, c[qp]) * params_.pressure * tr_strain_u * shape_test; 
                                    c_el_vec(j) +=  der_c_pres * dx(qp);
                                }

                                c_el_vec(j) += (elast * shape_test + frac) * dx(qp);
                                
                                if(params_.use_penalty_irreversibility){
                                    auto c_cold             = c[qp] - c_old[qp]; 
                                    auto c_cold_bracket     = c_cold < 0.0 ? c_cold: 0.0; 
                                    c_el_vec(j)             += params_.penalty_param * c_cold_bracket * shape_test * dx(qp);
                                }


                            }

                        }

                        U_view.add_vector(u_e, u_el_vec, g_view);
                        C_view.add_vector(c_e, c_el_vec, g_view);
                    }
                );
            }

            //check before boundary conditions
            if(check_derivatives_) {
                diff_ctrl_.check_grad(*this, x_const, g);
            }


            space_.apply_zero_constraints(g);

            // fully broken case is treated as Dirichlet BC
            if(params_.use_crack_set_irreversibiblity){
                apply_zero_constraints_irreversibiblity(g); 
            }

            // static int iter = 0;
            // write("g" + std::to_string(iter++) + ".m", g);
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override
        {
            if(empty(H)) {
                if(use_dense_hessian_) {
                    H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
                } else {
                    space_.create_matrix(H);
                }
            } else {
                H *= 0.0;
            }


            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

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

            //value based
            PrincipalStrains<USpace, Quadrature> strain(U, q);
            strain.update(x);

            //reference based
            PrincipalShapeStress<USpace, Quadrature> p_stress(U, q, params_.mu, params_.lambda);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

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

                //FIXME
                auto p_stress_view     = p_stress.view_device();

                auto H_view = space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
                        StaticMatrix<Scalar, Dim, Dim> strain_n, strain_p;
                        StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;
                        StaticMatrix<Scalar, Dim, Dim> stress; 

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

                        auto dx                 = differential_view.make(c_e);
                        auto c_grad_shape_el    = c_grad_shape_view.make(c_e);
                        auto c_shape_fun_el     = c_shape_view.make(c_e);

                        ////////////////////////////////////////////
                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {

                            const Scalar tr_strain_u = trace(el_strain.strain[qp]); 

                            const Scalar eep = elastic_energy(params_,  c[qp], tr_strain_u, el_strain.strain[qp]); 

                            for(SizeType l = 0; l < C_NDofs; ++l) {
                                const Scalar c_shape_l = c_shape_fun_el(l, qp);
                                auto &&      c_grad_l  = c_grad_shape_el(l, qp);

                                for(SizeType j = 0; j < C_NDofs; ++j) {

                                   el_mat(l, j) +=
                                        bilinear_cc(
                                            params_,
                                            c[qp],
                                            eep,
                                            c_shape_fun_el(j, qp),
                                            c_shape_l,
                                            c_grad_shape_el(j, qp),
                                            c_grad_l
                                    ) * dx(qp);        

                                    if(params_.pressure > 0){
                                        el_mat(l, j) +=  quadratic_degradation_deriv2(params_, c[qp]) * params_.pressure * tr_strain_u * c_shape_fun_el(j, qp) *  c_shape_l * dx(qp);        
                                    }

                                    if(params_.use_penalty_irreversibility){
                                        el_mat(l, j) += params_.penalty_param * c_shape_fun_el(j, qp) * c_shape_l * dx(qp);
                                    }                                        
                                }
                            }


                            for(SizeType l = 0; l < U_NDofs; ++l) {
                                auto &&u_grad_l = u_grad_shape_el(l, qp);

                                for(SizeType j = 0; j < U_NDofs; ++j) {
                                    el_mat(C_NDofs + l, C_NDofs + j) += bilinear_uu(
                                        params_,
                                        c[qp],
                                        p_stress_view.stress(j, qp), 
                                        u_grad_l
                                        ) * dx(qp);
                                }
                            }


                            //////////////////////////////////////////////////////////////////////////////////////////////////////
                            compute_stress(params_, trace(el_strain.strain[qp]), el_strain.strain[qp],  stress); 
                            // StaticMatrix<Scalar, Dim, Dim> identity; 
                            // identity.identity();                                         

                            for(SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                                //CHANGE (pre-compute/store shape fun)
                                const Scalar c_shape_i = c_shape_fun_el(c_i, qp);

                                for(SizeType u_i = 0; u_i < U_NDofs; ++u_i) {

                                    //CHANGE (use pre-computed strain)
                                    // auto strain_shape = 0.5 * (u_grad_shape_el(u_i, qp) + transpose(u_grad_shape_el(u_i, qp))); 
                                    auto && strain_shape = u_strain_shape_el(u_i, qp);

                                    Scalar val =
                                        bilinear_uc(
                                            params_,
                                             c[qp],
                                             stress,
                                             strain_shape,
                                             c_shape_i
                                        ) * dx(qp);

                                    if(params_.pressure > 0){
                                        const Scalar tr_strain_shape = sum(diag(strain_shape)); 
                                        val += quadratic_degradation_deriv(params_, c[qp]) * params_.pressure * tr_strain_shape * c_shape_i * dx(qp); 
                                    }

                                    el_mat(c_i, C_NDofs + u_i) += val;
                                    el_mat(C_NDofs + u_i, c_i) += val;
                                }
                            }



                        }

                        space_view.add_matrix(e, el_mat, H_view);
                    }
                );
            }

            //check before boundary conditions
            if(check_derivatives_) {
                diff_ctrl_.check_hessian(*this, x_const, H);
            }

            space_.apply_constraints(H);

            if(params_.use_crack_set_irreversibiblity){
                apply_zero_constraints_irreversibiblity(H);
            }

            // static int iter = 0;
            // write("H" + std::to_string(iter++) + ".m", H);
            return true;
        }

        //////////////////////////////////////////

        template<class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &elastic_energy,
            const Scalar &shape_trial,
            const Scalar &shape_test,
            const GradShape &grad_trial,
            const GradShape &grad_test
            )
        {
            return diffusion_c(params, grad_trial, grad_test) + reaction_c(params,  shape_trial,  shape_test) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy, shape_trial, shape_test);
        }

        // (sigma+(phi_u), epsilon(u)) * g'_c * phi_c
        template<class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cu(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Stress &stress_p,
            const FullStrain &full_strain,
            const Scalar &c_trial_fun
            )
        {
            return ((1.0 - params.regularization) * quadratic_degradation_deriv(params, phase_field_value)) * c_trial_fun * inner(stress_p, full_strain);
        }


        template<class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Stress &stress,
            const FullStrain &full_strain,
            const Scalar &c_trial_fun
            )
        {
            return quadratic_degradation_deriv(params, phase_field_value) * c_trial_fun * inner(stress, full_strain);
        }

        template<class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(
            const Parameters &params,
            const Scalar &phase_field_value,
            const StressShape &stress,
            const Grad &g_test)
        {
            const Scalar gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) + params.regularization);
            auto C_test  = 0.5 * (g_test  + transpose(g_test));
            return inner(gc * stress, C_test);
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
            const Scalar &elastic_energy,
            const Scalar &trial,
            const Scalar &test)
        {
            const Scalar dcc = (1.0 - params.regularization)*quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * trial * elastic_energy * test;
        }

        template<class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar grad_elastic_energy_wrt_c(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &trace,
            const Strain &strain)
        {
            return (quadratic_degradation_deriv(params, phase_field_value) * (1.0 - params.regularization))* strain_energy(params, trace, strain);
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

        template<class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(
            const Parameters &params,
            const Scalar &tr,
            const Strain &strain,
            Stress &stress)
        {
            stress = (2.0 * params.mu * strain) + (params.lambda * tr * (device::identity<Scalar>())); 
        }

    
        template<class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(
            const Parameters &params,
            //c
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad,
            // u
            const Scalar &trace,
            const Strain &strain)
        {
            return fracture_energy(params, phase_field_value, phase_field_grad) +
            elastic_energy(params,  phase_field_value, trace, strain);
        }

        template<class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar fracture_energy(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Grad   &phase_field_grad)
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
            return 0.5 * params.lambda * trace*trace + params.mu * inner(strain, strain);
        }

        template<class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(
            const Parameters &params,
            const Scalar &phase_field_value,
            const Scalar &trace,
            const Strain &strain)
        {
            return (quadratic_degradation(params, phase_field_value) * (1.0 - params.regularization) + params.regularization)* strain_energy(params, trace, strain);
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

        void old_solution(const Vector & x_old)
        {
            x_old_ = x_old; 
        }

        const Vector & old_solution() const 
        {
            return x_old_;
        }        

        void get_old_solution(Vector & x) const 
        {
            x = x_old_;
        }                


        void build_irreversility_constraint(Vector &lb)
        {   
            {
                auto d_x_old = const_device_view(x_old_);

                parallel_transform(lb, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                {
                    if(i%(Dim+1)==0){
                        return d_x_old.get(i); 
                    }
                    else{
                        return -9e15; 
                    }                    
                });
            }
        }    

        // this 2 functions need to be moved to BC conditions 
        void apply_zero_constraints_irreversibiblity(Vector & g) const
        {
            {
                auto d_x_old = const_device_view(x_old_);

                parallel_transform(g, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                {
                    if(i%(Dim+1)==0){
                        if(d_x_old.get(i) > params_.crack_set_tol){
                            return 0.0; 
                        }
                        else{
                            return xi; 
                        }
                    }
                    else{
                        return xi; 
                    }                    
                
                });
            }
        }

        // we should move this to BC conditions
        // also, will not run efficienetly in parallel 
        void apply_zero_constraints_irreversibiblity(Matrix & H) const
        {
            std::vector<SizeType> indices; 
            {
                Read<Vector> r(x_old_);

                Range range_w = range(x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(i%(Dim+1)==0 && x_old_.get(i) > params_.crack_set_tol){
                        indices.push_back(i);
                    }
                }
            }                 

            set_zero_rows(H, indices, 1.);
        }        


    private:
        FunctionSpace & space_;
        Parameters params_;
        DiffController<Matrix, Vector> diff_ctrl_;

        bool use_dense_hessian_;
        bool check_derivatives_;

        Vector x_old_; // stores old solution  - used for treatment of irreversibility constraint 
    };

}
#endif