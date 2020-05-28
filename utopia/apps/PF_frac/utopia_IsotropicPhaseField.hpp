#ifndef UTOPIA_PHASE_FIELD_HPP
#define UTOPIA_PHASE_FIELD_HPP

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

#include "utopia_petsc_NeumannBoundaryConditions.hpp"

#define UNROLL_FACTOR 4
#define U_MIN(a, b) ((a) < (b) ? (a) : (b))

namespace utopia {

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class IsotropicPhaseFieldForBrittleFractures final
        : public ExtendedFunction<typename FunctionSpace::Matrix, typename FunctionSpace::Vector> {
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

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        // using Quadrature = utopia::Quadrature<Shape, 2*(Shape::Order -1)>;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        class Parameters : public Configurable {
        public:
            void read(Input &in) override {
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
                in.get("use_pressure", use_pressure);

                in.get("use_penalty_irreversibility", use_penalty_irreversibility);
                in.get("penalty_param", penalty_param);

                in.get("use_crack_set_irreversibiblity", use_crack_set_irreversibiblity);
                in.get("crack_set_tol", crack_set_tol);

                in.get("mu", mu);
                in.get("lambda", lambda);
                in.get("fracture_toughness", fracture_toughness);
            }

            Parameters()
                : a(1.0),
                  b(1.0),
                  d(1.0),
                  f(1.0),
                  length_scale(0.0),
                  fracture_toughness(0.001),
                  mu(80.0),
                  lambda(120.0),
                  regularization(1e-10),
                  pressure(0.0),
                  penalty_param(0.0),
                  crack_set_tol(0.95)

            {}

            Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda;
            Scalar regularization, pressure, penalty_param, crack_set_tol;
            bool use_penalty_irreversibility{false}, use_crack_set_irreversibiblity{false}, use_pressure{false};
        };

        void read(Input &in) override {
            params_.read(in);
            in.get("use_dense_hessian", use_dense_hessian_);
            in.get("check_derivatives", check_derivatives_);
            in.get("diff_controller", diff_ctrl_);
            init_force_field(in);
        }

        // this is a bit of hack
        void set_pressure(const Scalar &pressure) { params_.pressure = pressure; }

        void use_crack_set_irreversibiblity(const bool &flg) { params_.use_crack_set_irreversibiblity = flg; }

        void init_force_field(Input &in) {
            in.get("neumann_bc", [&](Input &in) {
                in.get_all([&](Input &in) {
                    if (empty(force_field_)) {
                        space_.create_vector(force_field_);
                        force_field_.set(0.0);
                    }

                    NeumannBoundaryCondition<FunctionSpace> bc(space_);
                    bc.read(in);
                    bc.apply(force_field_);
                });
            });
        }

        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space)
            : space_(space), use_dense_hessian_(false), check_derivatives_(false) {
            if (params_.length_scale == 0) {
                params_.length_scale = 2.0 * space.mesh().min_spacing();
            }

            // this computation follows eq. 50 from "On penalization in variational phase-field models of britlle
            // fracture, Gerasimov, Lorenzis"
            if (params_.use_penalty_irreversibility) {
                Scalar tol = 1e-3;
                Scalar tol2 = tol * tol;
                params_.penalty_param = params_.fracture_toughness / params_.length_scale * (1.0 / tol2 - 1.0);
            }

            // in case of constant pressure field
            // if(params_.pressure){
            params_.use_pressure = true;
            setup_constant_pressure_field(params_.pressure);
            // }

            // needed for ML setup
            space_.create_vector(this->_x_eq_values);
            space_.create_vector(this->_eq_constrains_flg);

            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_pressure_field_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_pressure_field_);

            this->local_c_old_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_c_old_);
        }

        IsotropicPhaseFieldForBrittleFractures(FunctionSpace &space, const Parameters &params)
            : space_(space), params_(params), use_dense_hessian_(false), check_derivatives_(false) {
            this->local_x_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_x_);

            this->local_pressure_field_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_pressure_field_);

            this->local_c_old_ = std::make_shared<Vector>();
            space_.create_local_vector(*this->local_c_old_);
        }

        void use_dense_hessian(const bool val) { use_dense_hessian_ = val; }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
            // if(use_dense_hessian_) {
            //     H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
            // } else {
            space_.create_matrix(H);
            // }
            return true;
        }

        inline bool update(const Vector & /*x*/) override {
            // x_coeff_.update(x);
            return true;
        }

        bool value(const Vector &x_const, Scalar &val) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::value");

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            // auto &x = const_cast<Vector &>(x_const);

            // auto &x_old = const_cast<Vector &>(x_old_);
            // FEFunction<CSpace> c_old_fun(C, x_old);

            // auto &press = const_cast<Vector &>(pressure_field_);
            // FEFunction<CSpace> press_fun(C, press);

            // FEFunction<CSpace> c_fun(C, x);
            // FEFunction<USpace> u_fun(U, x);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            // udpate local pressure field
            space_.global_to_local(pressure_field_, *local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, local_pressure_field_);

            // update c_old
            space_.global_to_local(x_old_, *local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);
            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            val = 0.0;

            PrincipalStrains<USpace, Quadrature> strain(u_coeff, q);
            // PrincipalStrains<USpace, Quadrature> strain(U, q);
            // strain.update(x);

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

                Device::parallel_reduce(
                    space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        // StaticMatrix<Scalar, Dim, Dim> strain_n;
                        // StaticMatrix<Scalar, Dim, Dim> strain_p;

                        CElem c_e;
                        C_view.elem(i, c_e);

                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> c_old;
                        StaticVector<Scalar, NQuadPoints> p;
                        c_view.get(c_e, c);
                        c_old_view.get(c_e, c_old);
                        p_view.get(c_e, p);

                        UElem u_e;
                        U_view.elem(i, u_e);
                        auto el_strain = strain_view.make(u_e);
                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar el_energy = 0.0;

                        // #pragma clang loop unroll_count(U_MIN(NQuadPoints, UNROLL_FACTOR))
                        // #pragma GCC unroll U_MIN(NQuadPoints, UNROLL_FACTOR)
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            Scalar tr = trace(el_strain.strain[qp]);
                            if (params_.use_pressure) {
                                el_energy += quadratic_degradation(params_, c[qp]) * p[qp] * tr * dx(qp);
                            }

                            el_energy += energy(params_, c[qp], c_grad_el[qp], tr, el_strain.strain[qp]) * dx(qp);

                            if (params_.use_penalty_irreversibility) {
                                auto c_cold = c[qp] - c_old[qp];
                                auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                el_energy += params_.penalty_param / 2.0 * c_cold_bracket * c_cold_bracket * dx(qp);
                            }
                        }

                        assert(el_energy == el_energy);
                        return el_energy;
                    },
                    val);
            }

            val = x_const.comm().sum(val);

            assert(val == val);

            if (!empty(force_field_)) {
                // MAYBE -= dot(x_const, force_field_);
                val += dot(x_const, force_field_);
            }

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::value");
            return true;
        }

        bool gradient(const Vector &x_const, Vector &g) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::gradient");

            if (empty(g)) {
                space_.create_vector(g);
            } else {
                g.set(0.0);
            }

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            // udpate local pressure field
            space_.global_to_local(pressure_field_, *local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, local_pressure_field_);

            // update c_old
            space_.global_to_local(x_old_, *local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            // auto &x = const_cast<Vector &>(x_const);

            // auto &x_old = const_cast<Vector &>(x_old_);
            // FEFunction<CSpace> c_old_fun(C, x_old);

            // auto &press = const_cast<Vector &>(pressure_field_);
            // FEFunction<CSpace> press_fun(C, press);

            // FEFunction<CSpace> c_fun(C, x);
            // FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_old = c_old_fun.value(q);
            auto press_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            PrincipalStrains<USpace, Quadrature> strain(u_coeff, q);

            // PrincipalStrains<USpace, Quadrature> strain(U, q);
            // strain.update(x);

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

                auto g_view = space_.assembly_view_device(g);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    StaticMatrix<Scalar, Dim, Dim> stress;  //, strain_p;
                    StaticVector<Scalar, U_NDofs> u_el_vec;
                    StaticVector<Scalar, C_NDofs> c_el_vec;

                    u_el_vec.set(0.0);
                    c_el_vec.set(0.0);

                    ////////////////////////////////////////////

                    UElem u_e;
                    U_view.elem(i, u_e);
                    auto el_strain = strain_view.make(u_e);
                    // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
                    auto &&u_strain_shape_el = ref_strain_u_view.make(u_e);

                    ////////////////////////////////////////////

                    CElem c_e;
                    C_view.elem(i, c_e);
                    StaticVector<Scalar, NQuadPoints> c;
                    c_view.get(c_e, c);

                    StaticVector<Scalar, NQuadPoints> c_old;
                    c_old_view.get(c_e, c_old);

                    StaticVector<Scalar, NQuadPoints> p;
                    p_view.get(c_e, p);

                    auto c_grad_el = c_grad_view.make(c_e);
                    auto dx = differential_view.make(c_e);
                    auto c_grad_shape_el = c_grad_shape_view.make(c_e);
                    auto c_shape_fun_el = c_shape_view.make(c_e);

                    ////////////////////////////////////////////

                    // #pragma clang loop unroll_count(U_MIN(NQuadPoints, UNROLL_FACTOR))
                    // #pragma GCC unroll U_MIN(NQuadPoints, UNROLL_FACTOR)
                    for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                        const Scalar tr_strain_u = trace(el_strain.strain[qp]);

                        compute_stress(params_, tr_strain_u, el_strain.strain[qp], stress);
                        stress = (quadratic_degradation(params_, c[qp]) * (1.0 - params_.regularization) +
                                  params_.regularization) *
                                 stress;

                        // #pragma clang loop unroll_count(U_MIN(U_NDofs, UNROLL_FACTOR))
                        // #pragma GCC unroll U_MIN(U_NDofs, UNROLL_FACTOR)
                        for (SizeType j = 0; j < U_NDofs; ++j) {
                            auto &&strain_test = u_strain_shape_el(j, qp);
                            u_el_vec(j) += inner(stress, strain_test) * dx(qp);

                            if (params_.use_pressure) {
                                u_el_vec(j) +=
                                    quadratic_degradation(params_, c[qp]) * p[qp] * sum(diag(strain_test)) * dx(qp);
                            }
                        }

                        const Scalar elast =
                            grad_elastic_energy_wrt_c(params_, c[qp], tr_strain_u, el_strain.strain[qp]);

                        // #pragma clang loop unroll_count(U_MIN(C_NDofs, UNROLL_FACTOR))
                        // #pragma GCC unroll U_MIN(C_NDofs, UNROLL_FACTOR)
                        for (SizeType j = 0; j < C_NDofs; ++j) {
                            const Scalar shape_test = c_shape_fun_el(j, qp);
                            const Scalar frac = grad_fracture_energy_wrt_c(
                                params_, c[qp], c_grad_el[qp], shape_test, c_grad_shape_el(j, qp));

                            if (params_.use_pressure) {
                                const Scalar der_c_pres =
                                    quadratic_degradation_deriv(params_, c[qp]) * p[qp] * tr_strain_u * shape_test;
                                c_el_vec(j) += der_c_pres * dx(qp);
                            }

                            c_el_vec(j) += (elast * shape_test + frac) * dx(qp);

                            if (params_.use_penalty_irreversibility) {
                                auto c_cold = c[qp] - c_old[qp];
                                auto c_cold_bracket = c_cold < 0.0 ? c_cold : 0.0;
                                c_el_vec(j) += params_.penalty_param * c_cold_bracket * shape_test * dx(qp);
                            }
                        }
                    }

                    U_view.add_vector(u_e, u_el_vec, g_view);
                    C_view.add_vector(c_e, c_el_vec, g_view);
                });
            }

            // check before boundary conditions
            if (check_derivatives_) {
                diff_ctrl_.check_grad(*this, x_const, g);
            }

            if (!empty(force_field_)) {
                // MAYBE g -= force_field_;
                g += force_field_;
            }

            space_.apply_zero_constraints(g);

            // fully broken case is treated as Dirichlet BC
            if (params_.use_crack_set_irreversibiblity) {
                apply_zero_constraints_irreversibiblity(g, x_const);

                // // just a test...
                // auto* p_this = const_cast<IsotropicPhaseFieldForBrittleFractures<FunctionSpace> *>(this);
                // add_irr_values_markers(p_this->_x_eq_values, p_this->_eq_constrains_flg);
            }

            // static int iter = 0;
            // write("g" + std::to_string(iter++) + ".m", g);

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::gradient");
            return true;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override {
            UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::hessian");

            if (empty(H)) {
                // if(use_dense_hessian_) {
                //     H = local_zeros({space_.n_dofs(), space_.n_dofs()}); //FIXME
                // } else {
                space_.create_matrix(H);
                // }
            } else {
                H *= 0.0;
            }

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            // auto &x     = const_cast<Vector &>(x_const);t
            // auto &press = const_cast<Vector &>(pressure_field_);

            ////////////////////////////////////////////////////////////////////////////

            // update local vector x
            space_.global_to_local(x_const, *local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, local_x_);

            // udpate local pressure field
            space_.global_to_local(pressure_field_, *local_pressure_field_);
            auto p_coeff = std::make_shared<Coefficient<CSpace>>(C, local_pressure_field_);

            // update c_old
            space_.global_to_local(x_old_, *local_c_old_);
            auto c_old_coeff = std::make_shared<Coefficient<CSpace>>(C, local_c_old_);

            FEFunction<CSpace> c_old_fun(c_old_coeff);
            FEFunction<CSpace> press_fun(p_coeff);
            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto p_val = press_fun.value(q);

            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);
            auto differential = C.differential(q);

            // auto v_grad_shape = U.shape_grad(q);
            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            // value based
            PrincipalStrains<USpace, Quadrature> strain(u_coeff, q);
            // strain.update(x);

            // reference based
            PrincipalShapeStress<USpace, Quadrature> p_stress(U, q, params_.mu, params_.lambda);
            Strain<USpace, Quadrature> ref_strain_u(U, q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = space_.view_device();

                auto c_view = c_val.view_device();
                auto p_view = p_val.view_device();

                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                // auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                // FIXME
                auto p_stress_view = p_stress.view_device();

                auto H_view = space_.assembly_view_device(H);
                auto ref_strain_u_view = ref_strain_u.view_device();

                Device::parallel_for(space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    // StaticMatrix<Scalar, Dim, Dim> strain_n, strain_p;
                    StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;
                    StaticMatrix<Scalar, Dim, Dim> stress;

                    MixedElem e;
                    space_view.elem(i, e);
                    el_mat.set(0.0);

                    ////////////////////////////////////////////
                    UElem u_e;
                    U_view.elem(i, u_e);
                    auto el_strain = strain_view.make(u_e);
                    // auto u_grad_shape_el = v_grad_shape_view.make(u_e);
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
                    for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                        const Scalar tr_strain_u = trace(el_strain.strain[qp]);

                        const Scalar eep = elastic_energy(params_, c[qp], tr_strain_u, el_strain.strain[qp]);

                        // pragma GCCunroll(C_NDofs)
                        for (SizeType l = 0; l < C_NDofs; ++l) {
                            const Scalar c_shape_l = c_shape_fun_el(l, qp);
                            auto &&c_grad_l = c_grad_shape_el(l, qp);

                            // SYMMETRIC VERSION
                            for (SizeType j = l; j < C_NDofs; ++j) {
                                Scalar val = bilinear_cc(params_,
                                                         c[qp],
                                                         eep,
                                                         c_shape_fun_el(j, qp),
                                                         c_shape_l,
                                                         c_grad_shape_el(j, qp),
                                                         c_grad_l) *
                                             dx(qp);

                                if (params_.use_pressure) {
                                    val += quadratic_degradation_deriv2(params_, c[qp]) * p[qp] * tr_strain_u *
                                           c_shape_fun_el(j, qp) * c_shape_l * dx(qp);
                                }

                                if (params_.use_penalty_irreversibility) {
                                    val += params_.penalty_param * c_shape_fun_el(j, qp) * c_shape_l * dx(qp);
                                }

                                val = (l == j) ? (0.5 * val) : val;

                                el_mat(l, j) += val;
                                el_mat(j, l) += val;
                            }
                        }

                        // #pragma clang loop unroll_count(U_MIN(U_NDofs, UNROLL_FACTOR))
                        // #pragma GCC unroll U_MIN(U_NDofs, UNROLL_FACTOR)
                        for (SizeType l = 0; l < U_NDofs; ++l) {
                            auto &&u_strain_shape_l = u_strain_shape_el(l, qp);

                            // SYMMETRIC VERSION
                            el_mat(C_NDofs + l, C_NDofs + l) +=
                                bilinear_uu(params_, c[qp], p_stress_view.stress(l, qp), u_strain_shape_l) * dx(qp);

                            for (SizeType j = l + 1; j < U_NDofs; ++j) {
                                const Scalar v =
                                    bilinear_uu(params_, c[qp], p_stress_view.stress(j, qp), u_strain_shape_l) * dx(qp);

                                el_mat(C_NDofs + l, C_NDofs + j) += v;
                                el_mat(C_NDofs + j, C_NDofs + l) += v;
                            }
                        }

                        //////////////////////////////////////////////////////////////////////////////////////////////////////
                        compute_stress(params_, trace(el_strain.strain[qp]), el_strain.strain[qp], stress);

                        // #pragma clang loop unroll_count(U_MIN(C_NDofs, UNROLL_FACTOR))
                        // #pragma GCC unroll U_MIN(C_NDofs, UNROLL_FACTOR)
                        for (SizeType c_i = 0; c_i < C_NDofs; ++c_i) {
                            // CHANGE (pre-compute/store shape fun)
                            const Scalar c_shape_i = c_shape_fun_el(c_i, qp);

                            // #pragma clang loop unroll_count(U_MIN(U_NDofs, UNROLL_FACTOR))
                            // #pragma GCC unroll U_MIN(U_NDofs, UNROLL_FACTOR)
                            for (SizeType u_i = 0; u_i < U_NDofs; ++u_i) {
                                auto &&strain_shape = u_strain_shape_el(u_i, qp);

                                Scalar val = bilinear_uc(params_, c[qp], stress, strain_shape, c_shape_i) * dx(qp);

                                if (params_.use_pressure) {
                                    const Scalar tr_strain_shape = sum(diag(strain_shape));
                                    val += quadratic_degradation_deriv(params_, c[qp]) * p[qp] * tr_strain_shape *
                                           c_shape_i * dx(qp);
                                }

                                el_mat(c_i, C_NDofs + u_i) += val;
                                el_mat(C_NDofs + u_i, c_i) += val;
                            }
                        }
                    }

                    space_view.add_matrix(e, el_mat, H_view);
                });
            }

            // check before boundary conditions
            if (check_derivatives_) {
                diff_ctrl_.check_hessian(*this, x_const, H);
            }

            space_.apply_constraints(H);

            if (params_.use_crack_set_irreversibiblity) {
                apply_zero_constraints_irreversibiblity(H, x_const);
            }

            // static int iter = 0;
            // write("H" + std::to_string(iter++) + ".m", H);

            UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::hessian");
            return true;
        }

        //////////////////////////////////////////

        template <class GradShape>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Scalar &elastic_energy,
                                                         const Scalar &shape_trial,
                                                         const Scalar &shape_test,
                                                         const GradShape &grad_trial,
                                                         const GradShape &grad_test) {
            return diffusion_c(params, grad_trial, grad_test) + reaction_c(params, shape_trial, shape_test) +
                   elastic_deriv_cc(params, phase_field_value, elastic_energy, shape_trial, shape_test);
        }

        // (sigma+(phi_u), epsilon(u)) * g'_c * phi_c
        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_cu(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return ((1.0 - params.regularization) * quadratic_degradation_deriv(params, phase_field_value)) *
                   c_trial_fun * inner(stress_p, full_strain);
        }

        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return quadratic_degradation_deriv(params, phase_field_value) * c_trial_fun * inner(stress, full_strain);
        }

        // template<class StressShape, class Grad>
        // UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(
        //     const Parameters &params,
        //     const Scalar &phase_field_value,
        //     const StressShape &stress,
        //     const Grad &g_test)
        // {
        //     const Scalar gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) +
        //     params.regularization); auto C_test  = 0.5 * (g_test  + transpose(g_test)); return inner(gc * stress,
        //     C_test);
        // }

        template <class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(const Parameters &params,
                                                         const Scalar &phase_field_value,
                                                         const StressShape &stress,
                                                         const Grad &strain_test) {
            const Scalar gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) +
                               params.regularization);
            // auto C_test  = 0.5 * (g_test  + transpose(g_test));
            return inner(gc * stress, strain_test);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar diffusion_c(const Parameters &params,
                                                         const Grad &g_trial,
                                                         const Grad &g_test) {
            return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
        }

        UTOPIA_INLINE_FUNCTION static Scalar reaction_c(const Parameters &params,
                                                        const Scalar &trial,
                                                        const Scalar &test) {
            return (params.fracture_toughness / params.length_scale) * trial * test;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const Parameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &elastic_energy,
                                                              const Scalar &trial,
                                                              const Scalar &test) {
            const Scalar dcc = (1.0 - params.regularization) * quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * trial * elastic_energy * test;
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar grad_elastic_energy_wrt_c(const Parameters &params,
                                                                       const Scalar &phase_field_value,
                                                                       const Scalar &trace,
                                                                       const Strain &strain) {
            return (quadratic_degradation_deriv(params, phase_field_value) * (1.0 - params.regularization)) *
                   strain_energy(params, trace, strain);
        }

        template <class Grad, class GradTest>
        UTOPIA_INLINE_FUNCTION static Scalar grad_fracture_energy_wrt_c(const Parameters &params,
                                                                        const Scalar &phase_field_value,
                                                                        const Grad &phase_field_grad,
                                                                        const Scalar &test_function,
                                                                        const GradTest &grad_test_function) {
            return params.fracture_toughness * (1. / params.length_scale * phase_field_value * test_function +
                                                params.length_scale * inner(phase_field_grad, grad_test_function));
        }

        template <class Strain, class Stress>
        UTOPIA_INLINE_FUNCTION static void compute_stress(const Parameters &params,
                                                          const Scalar &tr,
                                                          const Strain &strain,
                                                          Stress &stress) {
            stress = (2.0 * params.mu * strain) + (params.lambda * tr * (device::identity<Scalar>()));
        }

        template <class Grad, class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar energy(const Parameters &params,
                                                    // c
                                                    const Scalar &phase_field_value,
                                                    const Grad &phase_field_grad,
                                                    // u
                                                    const Scalar &trace,
                                                    const Strain &strain) {
            return fracture_energy(params, phase_field_value, phase_field_grad) +
                   elastic_energy(params, phase_field_value, trace, strain);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar fracture_energy(const Parameters &params,
                                                             const Scalar &phase_field_value,
                                                             const Grad &phase_field_grad) {
            return params.fracture_toughness *
                   (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
                    params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar strain_energy(const Parameters &params,
                                                           const Scalar trace,
                                                           const Strain &strain) {
            return 0.5 * params.lambda * trace * trace + params.mu * inner(strain, strain);
        }

        template <class Strain>
        UTOPIA_INLINE_FUNCTION static Scalar elastic_energy(const Parameters &params,
                                                            const Scalar &phase_field_value,
                                                            const Scalar &trace,
                                                            const Strain &strain) {
            return (quadratic_degradation(params, phase_field_value) * (1.0 - params.regularization) +
                    params.regularization) *
                   strain_energy(params, trace, strain);
        }

        UTOPIA_INLINE_FUNCTION static Scalar degradation(const Parameters &params, const Scalar &c) {
            Scalar imc = 1.0 - c;
            Scalar res = params.f + imc * params.d;
            imc *= imc;
            res += params.b * imc;
            imc *= imc;  // FIXME

            res += params.a * imc;
            return res;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation(const Parameters &, const Scalar &c) {
            Scalar imc = 1.0 - c;
            return imc * imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv(const Parameters &, const Scalar &c) {
            Scalar imc = 1.0 - c;
            return -2.0 * imc;
        }

        UTOPIA_INLINE_FUNCTION static Scalar quadratic_degradation_deriv2(const Parameters &, const Scalar &) {
            return 2.0;
        }

        void old_solution(const Vector &x_old) { x_old_ = x_old; }

        Vector &old_solution() { return x_old_; }

        void get_old_solution(Vector &x) const { x = x_old_; }

        void set_old_solution(const Vector &x) { x_old_ = x; }

        void build_irreversility_constraint(Vector &lb) {
            {
                auto d_x_old = const_device_view(x_old_);

                auto lb_view = view_device(lb);
                parallel_for(range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                    if (i % (Dim + 1) == 0) {
                        lb_view.set(i, d_x_old.get(i));
                    } else {
                        lb_view.set(i, -9e15);
                    }
                });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void apply_zero_constraints_irreversibiblity(Vector &g) const {
            {
                auto d_x_old = const_device_view(x_old_);

                auto g_view = view_device(g);
                parallel_for(range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                    if (i % (Dim + 1) == 0) {
                        if (d_x_old.get(i) > params_.crack_set_tol) {
                            g_view.set(i, 0.0);
                        }
                    }
                });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const {
            {
                auto d_x_old = const_device_view(x_old_);
                auto d_x = const_device_view(x);

                auto g_view = view_device(g);
                parallel_for(range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                    if (i % (Dim + 1) == 0) {
                        if (d_x_old.get(i) > params_.crack_set_tol || d_x.get(i) > params_.crack_set_tol) {
                            g_view.set(i, 0.0);
                        }
                    }
                });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void add_irr_values_markers(Vector &val, Vector &flg) const {
            {
                auto d_x_old = const_device_view(x_old_);

                auto val_view = view_device(val);
                parallel_for(range_device(val), UTOPIA_LAMBDA(const SizeType &i) {
                    if (i % (Dim + 1) == 0) {
                        if (d_x_old.get(i) > params_.crack_set_tol) {
                            val_view.set(i, d_x_old.get(i));
                        }
                    }
                });

                auto flg_view = view_device(flg);
                parallel_for(range_device(flg), UTOPIA_LAMBDA(const SizeType &i) {
                    if (i % (Dim + 1) == 0) {
                        if (d_x_old.get(i) > params_.crack_set_tol) {
                            flg_view.set(i, 1.0);
                        }
                    }
                });
            }
        }

        // we should move this to BC conditions
        // also, will not run efficienetly in parallel
        void apply_zero_constraints_irreversibiblity(Matrix &H) const {
            std::vector<SizeType> indices;
            {
                Read<Vector> r(x_old_);

                Range range_w = range(x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 0 && x_old_.get(i) > params_.crack_set_tol) {
                        indices.push_back(i);
                    }
                }
            }

            set_zero_rows(H, indices, 1.);
        }

        void apply_zero_constraints_irreversibiblity(Matrix &H, const Vector &x) const {
            std::vector<SizeType> indices;
            {
                Read<Vector> r(x_old_);
                Read<Vector> r2(x);

                Range range_w = range(x_old_);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 0 &&
                        (x_old_.get(i) > params_.crack_set_tol || x.get(i) > params_.crack_set_tol)) {
                        indices.push_back(i);
                    }
                }
            }

            set_zero_rows(H, indices, 1.);
        }

        void pressure_field(const Vector &pressure_field) { pressure_field_ = pressure_field; }

        void setup_constant_pressure_field(const Scalar &p_val) {
            if (empty(pressure_field_)) {
                space_.create_vector(pressure_field_);
            }

            pressure_field_.set(p_val);
        }

        Vector &pressure_field() {
            if (empty(pressure_field_)) space_.create_vector(pressure_field_);

            return pressure_field_;
        }

    private:
        FunctionSpace &space_;
        Parameters params_;
        DiffController<Matrix, Vector> diff_ctrl_;

        bool use_dense_hessian_;
        bool check_derivatives_;

        Vector x_old_;           // stores old solution  - used for treatment of irreversibility constraint
        Vector pressure_field_;  // stores heterogenous pressure field - ideally, this vector would have lower size than
                                 // all 4 variables

        Vector force_field_;

        std::shared_ptr<Vector> local_x_;
        std::shared_ptr<Vector> local_pressure_field_;
        std::shared_ptr<Vector> local_c_old_;
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
