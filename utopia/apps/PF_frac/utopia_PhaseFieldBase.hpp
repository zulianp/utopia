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

    template <class FunctionSpace>
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

            in.get("turn_off_uc_coupling", turn_off_uc_coupling);
            in.get("turn_off_cu_coupling", turn_off_cu_coupling);
        }

        PFFracParameters()
            : a(1.0),
              b(1.0),
              d(1.0),
              f(1.0),
              length_scale(0.0),
              fracture_toughness(0.001),
              mu(80.0),
              lambda(120.0),
              regularization(1e-15),
              pressure(0.0),
              penalty_param(0.0),
              crack_set_tol(0.93)

        {}

        Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda;
        Scalar regularization, pressure, penalty_param, crack_set_tol;
        bool use_penalty_irreversibility{false}, use_crack_set_irreversibiblity{false}, use_pressure{false};
        bool turn_off_uc_coupling{false}, turn_off_cu_coupling{false};
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
        // using Quadrature = utopia::Quadrature<Shape, 2*(Shape::Order -1)>;
        using Quadrature = utopia::Quadrature<Shape, 2 * (Shape::Order)>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        void read(Input &in) override {
            params_.read(in);
            // in.get("use_dense_hessian", use_dense_hessian_);
            // in.get("check_derivatives", check_derivatives_);
            // in.get("diff_controller", diff_ctrl_);
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

        PhaseFieldFracBase(FunctionSpace &space) : space_(space) {
            if (params_.length_scale == 0) {
                params_.length_scale = 2.0 * space.mesh().min_spacing();
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

        //////////////////////////////////////////

        // template <class GradShape>
        // UTOPIA_INLINE_FUNCTION static Scalar bilinear_cc(const PFFracParameters &params,
        //                                                  const Scalar &phase_field_value,
        //                                                  const Scalar &elastic_energy,
        //                                                  const Scalar &shape_trial,
        //                                                  const Scalar &shape_test,
        //                                                  const GradShape &grad_trial,
        //                                                  const GradShape &grad_test) {
        //     return diffusion_c(params, grad_trial, grad_test) + reaction_c(params, shape_trial, shape_test) +
        //            elastic_deriv_cc(params, phase_field_value, elastic_energy, shape_trial, shape_test);
        // }

        template <class StressShape, class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar bilinear_uu(const PFFracParameters &params,
                                                         const Scalar &phase_field_value,
                                                         const StressShape &stress,
                                                         const Grad &strain_test) {
            const Scalar gc = ((1.0 - params.regularization) * quadratic_degradation(params, phase_field_value) +
                               params.regularization);
            return inner(gc * stress, strain_test);
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

        // template <class Grad, class Strain>
        // UTOPIA_INLINE_FUNCTION static Scalar energy(const PFFracParameters &params,
        //                                             // c
        //                                             const Scalar &phase_field_value,
        //                                             const Grad &phase_field_grad,
        //                                             // u
        //                                             const Scalar &trace,
        //                                             const Strain &strain) {
        //     return fracture_energy(params, phase_field_value, phase_field_grad) +
        //            elastic_energy(params, phase_field_value, trace, strain);
        // }

        // template <class Grad>
        // UTOPIA_INLINE_FUNCTION static Scalar fracture_energy(const PFFracParameters &params,
        //                                                      const Scalar &phase_field_value,
        //                                                      const Grad &phase_field_grad) {
        //     return params.fracture_toughness *
        //            (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
        //             params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
        // }

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

        void old_solution(const Vector &x_old) { x_old_ = x_old; }

        Vector &old_solution() { return x_old_; }

        void get_old_solution(Vector &x) const { x = x_old_; }

        void set_old_solution(const Vector &x) { x_old_ = x; }

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
        void apply_zero_constraints_irreversibiblity(Vector &g, const Vector &x) const {
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

    protected:
        FunctionSpace &space_;
        PFFracParameters params_;
        DiffController<Matrix, Vector> diff_ctrl_;

        bool use_dense_hessian_ = {false};
        bool check_derivatives_ = {false};

        Vector x_old_;           // stores old solution  - used for treatment of
                                 // irreversibility constraint
        Vector pressure_field_;  // stores heterogenous pressure field - ideally, this
                                 // vector would have lower size than all 4 variables

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
