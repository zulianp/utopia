#ifndef UTOPIA_PHASE_FIELD_BASE_HPP
#define UTOPIA_PHASE_FIELD_BASE_HPP

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
    class PFFracParameters : public Configurable {
    public:
        using Scalar = typename FunctionSpace::Scalar;
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
            in.get("nu", nu);
            in.get("E", E);
            in.get("l_0", l_0);
            in.get("pressure0", pressure0);

            in.get("turn_off_uc_coupling", turn_off_uc_coupling);
            in.get("turn_off_cu_coupling", turn_off_cu_coupling);

            in.get("mobility", mobility);
            in.get("use_mobility", use_mobility);

            kappa = lambda + (2.0 * mu / Dim);

            if (nu != 0.0 && E != 0.0) {
                mu = E / (2.0 * (1. + nu));
                lambda = (2.0 * nu * mu) / (1.0 - (2.0 * nu));

                if (mpi_world_rank() == 0) {
                    utopia::out() << "mu: " << mu << "  lambda: " << lambda << "  Gc: " << fracture_toughness << "  \n";
                }
            }
        }

        PFFracParameters()
            : a(1.0),
              b(1.0),
              d(1.0),
              f(1.0),
              length_scale(0.0),
              fracture_toughness(1e-3),
              mu(80.0),
              lambda(120.0),
              // mu(100.0),
              // lambda(100.0),
              nu(0.0),
              E(0.0),
              l_0(1.0),
              pressure0(1e-3),
              regularization(1e-10),
              pressure(0.0),
              penalty_param(0.0),
              crack_set_tol(0.93),
              // mobility(1e-5)
              mobility(1e-6)

        {
            kappa = lambda + (2.0 * mu / Dim);
        }

        bool kroneckerDelta(const SizeType &i, const SizeType &j) { return (i == j) ? 1.0 : 0.0; }

        void fill_in_isotropic_elast_tensor() {
            for (SizeType i = 0; i < Dim; ++i) {
                for (SizeType j = 0; j < Dim; ++j) {
                    for (SizeType k = 0; k < Dim; ++k) {
                        for (SizeType l = 0; l < Dim; ++l) {
                            Scalar val = this->lambda * kroneckerDelta(i, j) * kroneckerDelta(k, l);
                            val += this->mu * (kroneckerDelta(i, k) * kroneckerDelta(j, l));
                            val += this->mu * (kroneckerDelta(i, l) * kroneckerDelta(j, k));
                            elast_tensor.set(i, j, k, l, val);
                        }
                    }
                }
            }

            I4sym.identity_sym();
            kappa = lambda + (2.0 * mu / Dim);
        }

        Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda, kappa, nu, E, l_0, pressure0;
        Scalar regularization, pressure, penalty_param, crack_set_tol, mobility;
        bool use_penalty_irreversibility{false}, use_crack_set_irreversibiblity{false}, use_pressure{false};
        bool turn_off_uc_coupling{false}, turn_off_cu_coupling{false};
        bool use_mobility{false};

        Tensor4th<Scalar, Dim, Dim, Dim, Dim> elast_tensor;
        Tensor4th<Scalar, Dim, Dim, Dim, Dim> I4sym;
    };

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class PhaseFieldFracBase : public ExtendedFunction<typename FunctionSpace::Matrix, typename FunctionSpace::Vector> {
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

        using PFFracParameters = utopia::PFFracParameters<FunctionSpace>;

        // FIXME
        using Shape = typename FunctionSpace::Shape;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        // using ExtendedFunction<typename FunctionSpace::Matrix, typename
        // FunctionSpace::Vector>::get_eq_constrains_flg; using
        // ExtendedFunction<typename FunctionSpace::Matrix,
        //                        typename
        //                        FunctionSpace::Vector>::get_eq_constrains_values;

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
        void use_penalty_irreversibility(const bool &flg) { params_.use_penalty_irreversibility = flg; }

        void turn_off_uc_coupling(const bool &flg) { params_.turn_off_uc_coupling = flg; }
        void turn_off_cu_coupling(const bool &flg) { params_.turn_off_cu_coupling = flg; }

        void init_force_field(Input &in) {
            in.get("neumann_bc", [&](Input &in) {
                in.get_all([&](Input & /*in*/) {
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

        PhaseFieldFracBase(FunctionSpace &space) : space_(space) {
            if (params_.length_scale == 0) {
                params_.length_scale = 2.0 * space.mesh().min_spacing();
                // params_.length_scale = 2.0;
                if (mpi_world_rank() == 0) {
                    utopia::out() << "using ls = " << params_.length_scale << "  \n";
                }
            }

            // this computation follows eq. 50 from "On penalization in variational
            // phase-field models of britlle fracture, Gerasimov, Lorenzis"
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

        PhaseFieldFracBase(FunctionSpace &space, const PFFracParameters &params) : space_(space), params_(params) {
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

        ////////////////////////////////////////////////////////////////////////////////////
        virtual bool fracture_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;
        virtual bool elastic_energy(const Vector & /*x_const*/, Scalar & /*val*/) const = 0;

        // compute total crack volume (TCV), so we can compare to exact solution
        virtual bool compute_tcv(const Vector &x_const, Scalar &error) const {
            UTOPIA_TRACE_REGION_BEGIN("PFBase::compute_tcv");
            const Scalar PI = 3.141592653589793238463;

            Scalar tcv_exact = 0.0;
            Scalar computed_tcv = 0.0;

            USpace U;
            this->space_.subspace(1, U);
            CSpace C = this->space_.subspace(0);

            CSpace U1 = this->space_.subspace(1);
            CSpace U2 = this->space_.subspace(2);
            CSpace U3 = this->space_.subspace(3);

            ///////////////////////////////////////////////////////////////////////////

            // update local vector x
            this->space_.global_to_local(x_const, *this->local_x_);
            auto u_coeff = std::make_shared<Coefficient<USpace>>(U, this->local_x_);
            auto c_coeff = std::make_shared<Coefficient<CSpace>>(C, this->local_x_);

            auto u_coeff1 = std::make_shared<Coefficient<CSpace>>(U1, this->local_x_);
            auto u_coeff2 = std::make_shared<Coefficient<CSpace>>(U2, this->local_x_);
            auto u_coeff3 = std::make_shared<Coefficient<CSpace>>(U3, this->local_x_);

            FEFunction<CSpace> c_fun(c_coeff);
            FEFunction<USpace> u_fun(u_coeff);

            FEFunction<CSpace> u1_fun(u_coeff1);
            FEFunction<CSpace> u2_fun(u_coeff2);
            FEFunction<CSpace> u3_fun(u_coeff3);

            ////////////////////////////////////////////////////////////////////////////

            Quadrature q;

            auto c_val = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val = u_fun.value(q);

            auto u1_val = u1_fun.value(q);
            auto u2_val = u2_fun.value(q);
            auto u3_val = u3_fun.value(q);

            auto differential = C.differential(q);

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto U1_view = U1.view_device();
                auto U2_view = U2.view_device();
                auto U3_view = U3.view_device();

                auto c_view = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view = u_val.view_device();

                auto u1_view = u1_val.view_device();
                auto u2_view = u2_val.view_device();
                auto u3_view = u3_val.view_device();

                auto differential_view = differential.view_device();

                Device::parallel_reduce(
                    this->space_.element_range(),
                    UTOPIA_LAMBDA(const SizeType &i) {
                        StaticVector<Scalar, NQuadPoints> c;
                        StaticVector<Scalar, NQuadPoints> u1;
                        StaticVector<Scalar, NQuadPoints> u2;
                        StaticVector<Scalar, NQuadPoints> u3;

                        CElem c_e;
                        C_view.elem(i, c_e);
                        c_view.get(c_e, c);

                        CElem u1_e;
                        U1_view.elem(i, u1_e);
                        u1_view.get(u1_e, u1);

                        CElem u2_e;
                        U2_view.elem(i, u2_e);
                        u2_view.get(u2_e, u2);

                        CElem u3_e;
                        U3_view.elem(i, u3_e);
                        u3_view.get(u3_e, u3);

                        auto c_grad_el = c_grad_view.make(c_e);

                        auto dx = differential_view.make(c_e);

                        Scalar tcv_mine = 0.0;

                        // TCV = \int_\Omega u \cdot \nabla \phi
                        for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            if (Dim == 2) {
                                tcv_mine += ((u1[qp] * c_grad_el[qp](0)) + (u2[qp] * c_grad_el[qp](1))) * dx(qp);
                            } else {
                                tcv_mine += ((u1[qp] * c_grad_el[qp](0)) + (u2[qp] * c_grad_el[qp](1)) +
                                             (u3[qp] * c_grad_el[qp](2))) *
                                            dx(qp);
                            }
                        }

                        assert(tcv_mine == tcv_mine);
                        return tcv_mine;
                    },
                    computed_tcv);
            }

            computed_tcv = x_const.comm().sum(computed_tcv);
            assert(computed_tcv == computed_tcv);

            // params depends on initial fracture lenght and pressure
            if (Dim == 2) {
                tcv_exact = 2.0 * this->params_.pressure0 * this->params_.l_0 * this->params_.l_0 *
                            (1.0 - this->params_.nu * this->params_.nu) * PI / this->params_.E;
            } else {
                tcv_exact = 16.0 * this->params_.pressure0 * this->params_.l_0 * this->params_.l_0 * this->params_.l_0 *
                            (1.0 - this->params_.nu * this->params_.nu) / this->params_.E / 3.0;
            }

            error = device::abs(computed_tcv - tcv_exact);
            if (mpi_world_rank() == 0) {
                std::cout << "computed_tcv: " << computed_tcv << "  exact: " << tcv_exact << "  error: " << error
                          << "\n ";
            }

            UTOPIA_TRACE_REGION_END("PFBase::compute_tcv");
            return true;
        }

        virtual void compute_cod(const Vector &x_const, Scalar &error) const {
            const Scalar PI = 3.141592653589793238463;

            // coordinates of the point at which we compute crack openning
            Scalar x_coord = 0.0;
            Scalar y_coord = 0.0;

            Scalar cod_exact = 4.0 / PI * this->params_.pressure0 * this->params_.l_0 *
                               (1. - (this->params_.nu * this->params_.nu)) / this->params_.E;

            Scalar rho = std::sqrt((x_coord * x_coord) + (y_coord * y_coord));
            Scalar rho_l0 = (rho / this->params_.l_0);

            cod_exact = cod_exact * std::sqrt(1. - rho_l0 * rho_l0);

            // search for max disp_y (should be symmetric)
            Scalar cod_computed = 0.0;

            {
                Read<Vector> r(x_const);

                Range range_w = range(x_const);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++) {
                    if (i % (Dim + 1) == 2 && x_const.get(i) > cod_computed) {
                        cod_computed = x_const.get(i);
                    }
                }
            }

            cod_computed = x_const.comm().max(cod_computed);

            error = device::abs(cod_computed - cod_exact);

            if (mpi_world_rank() == 0) {
                std::cout << "cod_exact: " << cod_exact << "  cod_computed: " << cod_computed << "  error: " << error
                          << "\n ";
            }
        }

        virtual void update_history_field(const Vector & /*x_const*/) const {}

        ////////////////////////////////////////////////////////////////////////////////////

        template <typename PhaseFieldValue, class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uu(const PFFracParameters &params,
                                                                  const PhaseFieldValue &phase_field_value,
                                                                  const StressShape &stress,
                                                                  const Grad &strain_test) {
            const auto gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) +
                             params.regularization);
            return inner(gc * stress, strain_test);
        }

        template <class Grad>
        UTOPIA_INLINE_FUNCTION static auto diffusion_c(const PFFracParameters &params,
                                                       const Grad &g_trial,
                                                       const Grad &g_test) {
            return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
        }

        template <typename FunValue>
        UTOPIA_INLINE_FUNCTION static FunValue reaction_c(const PFFracParameters &params,
                                                          const FunValue &trial,
                                                          const FunValue &test) {
            return (params.fracture_toughness / params.length_scale) * trial * test;
        }

        UTOPIA_INLINE_FUNCTION static Scalar elastic_deriv_cc(const PFFracParameters &params,
                                                              const Scalar &phase_field_value,
                                                              const Scalar &elastic_energy_positive,
                                                              const Scalar &trial,
                                                              const Scalar &test) {
            const Scalar dcc = (1.0 - params.regularization) * quadratic_degradation_deriv2(params, phase_field_value);
            return dcc * trial * elastic_energy_positive * test;
        }

        template <typename PhaseFieldValue, class Grad, typename TestFunction, class GradTest>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue grad_fracture_energy_wrt_c(
            const PFFracParameters &params,
            const PhaseFieldValue &phase_field_value,
            const Grad &phase_field_grad,
            const TestFunction &test_function,
            const GradTest &grad_test_function) {
            return params.fracture_toughness * ((1. / params.length_scale * phase_field_value * test_function) +
                                                (params.length_scale * inner(phase_field_grad, grad_test_function)));
        }

        template <class Stress, class FullStrain>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uc(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const Stress &stress_p,
                                                         const FullStrain &full_strain,
                                                         const Scalar &c_trial_fun) {
            return c_trial_fun * inner(quadratic_degradation_deriv(params, phase_field_value) * stress_p, full_strain);
        }

        template <typename PhaseFieldValue, class Grad>
        UTOPIA_INLINE_FUNCTION static PhaseFieldValue fracture_energy(const PFFracParameters &params,
                                                                      const PhaseFieldValue &phase_field_value,
                                                                      const Grad &phase_field_grad) {
            return params.fracture_toughness *
                   (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
                    params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
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

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation(const PFFracParameters &, const C &c) {
            C imc = 1.0 - c;
            return imc * imc;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation_deriv(const PFFracParameters &, const C &c) {
            C imc = 1.0 - c;
            return -2.0 * imc;
        }

        template <typename C>
        UTOPIA_INLINE_FUNCTION static C quadratic_degradation_deriv2(const PFFracParameters &, const C &) {
            return 2.0;
        }

        Vector &old_solution() { return x_old_; }

        void get_old_solution(Vector &x) const { x = x_old_; }

        void set_old_solution(const Vector &x) {
            x_old_ = x;
            update_history_field(x_old_);
        }

        void build_irreversility_constraint(Vector &lb) {
            {
                auto d_x_old = const_device_view(x_old_);

                auto lb_view = view_device(lb);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            lb_view.set(i, d_x_old.get(i));
                        } else {
                            lb_view.set(i, -9e15);
                        }
                    });
            }
        }

        void build_irreversility_constraint(Vector &lb, Vector &ub) {
            {
                auto d_x_old = const_device_view(x_old_);

                auto lb_view = view_device(lb);
                auto ub_view = view_device(ub);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            lb_view.set(i, d_x_old.get(i));
                            ub_view.set(i, 1.0);
                        } else {
                            lb_view.set(i, -9e15);
                            ub_view.set(i, 9e15);
                        }
                    });
            }
        }

        void make_iterate_feasible(const Vector &lb, const Vector &ub, Vector &x) {
            {
                auto d_x_old = view_device(x);

                auto lb_view = const_device_view(lb);
                auto ub_view = const_device_view(ub);
                parallel_for(
                    range_device(lb), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            // Scalar li = lb_view.get(i);
                            Scalar ui = ub_view.get(i);
                            auto xi = d_x_old.get(i);
                            // if (li >= xi) {
                            //     d_x_old.set(i, li);
                            // }
                            if (ui <= xi) {
                                d_x_old.set(i, ui);
                            }
                            // else {
                            //     // d_x_old.set(i, (ui <= xi) ? ui : xi);
                            // }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void apply_zero_constraints_irreversibiblity(Vector &g) const {
            {
                auto d_x_old = const_device_view(x_old_);

                auto g_view = view_device(g);
                parallel_for(
                    range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                g_view.set(i, 0.0);
                            }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        virtual void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const {
            {
                auto d_x_old = const_device_view(x_old_);
                auto d_x = const_device_view(x);

                auto g_view = view_device(g);
                parallel_for(
                    range_device(g), UTOPIA_LAMBDA(const SizeType &i) {
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
                parallel_for(
                    range_device(val), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                val_view.set(i, d_x_old.get(i));
                            }
                        }
                    });

                auto flg_view = view_device(flg);
                parallel_for(
                    range_device(flg), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                flg_view.set(i, 1.0);
                            }
                        }
                    });
            }
        }

        // this 2 functions need to be moved to BC conditions
        void add_irr_values_markers(Vector &val, Vector &flg, const Vector &x_current) const {
            {
                auto d_x_old = const_device_view(x_current);

                auto val_view = view_device(val);
                parallel_for(
                    range_device(val), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                val_view.set(i, d_x_old.get(i));
                            }
                        }
                    });

                auto flg_view = view_device(flg);
                parallel_for(
                    range_device(flg), UTOPIA_LAMBDA(const SizeType &i) {
                        if (i % (Dim + 1) == 0) {
                            if (d_x_old.get(i) > params_.crack_set_tol) {
                                flg_view.set(i, 1.0);
                            }
                        }
                    });
            }
        }

        void add_pf_constraints(const Vector &x_current) const {
            auto *p_this = const_cast<PhaseFieldFracBase<FunctionSpace> *>(this);

            Vector &bc_flgs = p_this->get_eq_constrains_flg();
            Vector &bc_values = p_this->get_eq_constrains_values();
            p_this->space_.apply_constraints(bc_values);
            p_this->space_.build_constraints_markers(bc_flgs);
            p_this->add_irr_values_markers(bc_values, bc_flgs, x_current);
        }

        // we should move this to BC conditions
        // also, will not run efficienetly in parallel
        virtual void apply_zero_constraints_irreversibiblity(Matrix &H) const {
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

        virtual void apply_zero_constraints_irreversibiblity(Matrix &H, const Vector &x) const {
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

        void set_dt(const Scalar &dt) { dt_ = dt; }

        Scalar get_dt() const { return dt_; }

        virtual void write_to_file(const std::string &output_path, const Vector &x, const Scalar time) {
            space_.write(output_path + "_" + std::to_string(time) + ".vtr", x);
        }

        virtual bool must_reduce_time_step(const Vector &) {
            return false;
        }

    protected:
        FunctionSpace &space_;
        PFFracParameters params_;
        DiffController<Matrix, Vector> diff_ctrl_;

        bool use_dense_hessian_{false};
        bool check_derivatives_{false};

        Vector x_old_;           // stores old solution  - used for treatment of
                                 // irreversibility constraint
        Vector pressure_field_;  // stores heterogenous pressure field - ideally, this
                                 // vector would have lower size than all 4 variables

        Vector force_field_;

        std::shared_ptr<Vector> local_x_;
        std::shared_ptr<Vector> local_pressure_field_;
        std::shared_ptr<Vector> local_c_old_;

        Scalar dt_;
    };

}  // namespace utopia

// clean-up macros
#undef UNROLL_FACTOR
#undef U_MIN
#endif
