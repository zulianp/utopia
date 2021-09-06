#ifndef UTOPIA_KOKKOS_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP
#define UTOPIA_KOKKOS_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP

#include "utopia_Tracer.hpp"
#include "utopia_Views.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_PhaseFieldKernels.hpp"
#include "utopia_kokkos_Strain.hpp"
#include "utopia_kokkos_SubView.hpp"

#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_kokkos_StressLinearElasticityOp.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, int Dim>
        class IsotropicPhaseFieldForBrittleFractures : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using Gradient = typename FE::Gradient;
            using Function = typename FE::Function;
            using Measure = typename FE::Measure;
            using DynRankView = typename FE::DynRankView;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using VectorView = typename Super::VectorView;

            static constexpr int PHASE_FIELD_OFFSET = 0;
            static constexpr int DISPLACEMENT_OFFSET = 1;

            class Params : public Configurable {
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
                    }

                    in.get("displacement", displacement);
                    in.get("phase_field", phase_field);
                }

                Params()
                    : a(1.0),
                      b(1.0),
                      d(1.0),
                      f(1.0),
                      length_scale(0.0),
                      fracture_toughness(1e-3),
                      mu(80.0),
                      lambda(120.0),
                      nu(0.0),
                      E(0.0),
                      l_0(1.0),
                      pressure0(1e-3),
                      regularization(1e-10),
                      pressure(0.0),
                      penalty_param(0.0),
                      crack_set_tol(0.93),
                      mobility(1e-6) {
                    kappa = lambda + (2.0 * mu / Dim);
                }

                Scalar a, b, d, f, length_scale, fracture_toughness, mu, lambda, kappa, nu, E, l_0, pressure0;
                Scalar regularization, pressure, penalty_param, crack_set_tol, mobility;
                bool use_penalty_irreversibility{false};
                bool use_crack_set_irreversibiblity{false};
                bool use_pressure{false};
                bool turn_off_uc_coupling{false};
                bool turn_off_cu_coupling{false};
                bool use_mobility{false};
                int displacement{0};
                int phase_field{Dim};
            };

            using QuadraticDegradation = utopia::kokkos::QuadraticDegradation<Params>;
            using QuadraticPhaseFieldKernels = utopia::kokkos::PhaseFieldKernels<Params, QuadraticDegradation>;

            using V = StaticVector<Scalar, Dim>;
            using M = StaticMatrix<Scalar, Dim, Dim>;
            using SIMDType = Scalar;

            using InterpolateField = typename utopia::kokkos::Field<FE>::Interpolate;
            using InterpolateStrain =
                typename utopia::kokkos::kernels::InterpolateLinearizedStrainOp<Dim, Scalar, Gradient, DynRankView>;
            using InterpolatePhaseFieldGradient = typename utopia::kokkos::Gradient<FE>::Rank1Op;

            using GradUView = utopia::kokkos::SubView<DynRankView,
                                                      IdentityRange,
                                                      IdentityRange,
                                                      IdentityRange,
                                                      StaticRange<DISPLACEMENT_OFFSET, DISPLACEMENT_OFFSET + Dim>>;

            using LinearElasticityOp = utopia::kokkos::kernels::
                LinearElasticityOp<Dim, Scalar, Scalar, Scalar, Gradient, typename FE::Measure>;

            using CoeffLinearElasticityOp = utopia::kokkos::kernels::
                CoeffLinearElasticityOp<Dim, Scalar, Scalar, Scalar, DynRankView, typename FE::Measure>;

            IsotropicPhaseFieldForBrittleFractures(const std::shared_ptr<FE> &fe, Params op = Params())
                : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim + 1; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return true; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "IsotropicPhaseFieldForBrittleFractures"; }

            virtual bool update(const std::shared_ptr<Field<FE>> &x) override {
                if (!Super::update(x)) {
                    return false;
                }

                assert(x);
                assert(x->is_coefficient());

                if (!x->is_coefficient()) {
                    Utopia::Abort("IsotropicPhaseFieldForBrittleFractures::update, x must me in coefficient form!");
                }

                if (!gradient_) {
                    // Initialize gradient
                    gradient_ = std::make_shared<utopia::kokkos::Gradient<FE>>(this->fe_ptr());
                }

                gradient_->init(*x);
                return true;
            }

            class VectorOp {};

            class MatrixOp {};

            bool assemble_scalar() override {
                using DegradationFunction = utopia::kokkos::QuadraticDegradationKernel<InterpolateField, Scalar>;

                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_scalar");

                this->ensure_scalar_accumulator();
                auto data = this->scalar_data();

                // FE
                auto &fe = this->fe();
                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();
                const int n_quad_points = fe.n_quad_points();

                // Fields
                const int displacement = op_.displacement;
                const int phase_field = op_.phase_field;

                auto x = this->current_solution()->interpolate();
                auto p = pressure_field_->interpolate();
                auto grad_x = gradient_->data();

                // Model parameters
                const Scalar regularization = op_.regularization;
                const Scalar fracture_toughness = op_.fracture_toughness;
                const Scalar length_scale = op_.length_scale;

                // Function and operators
                DegradationFunction degradation_function(x, phase_field);
                CoeffLinearElasticityOp stress_op(op_.lambda, op_.mu, grad_x, measure, displacement);

                this->loop_cell(
                    "IsotropicPhaseFieldForBrittleFractures::assemble_scalar", UTOPIA_LAMBDA(const int &cell) {
                        Scalar integr = 0;
                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            Scalar val = 0.0;

                            //////////////////////////////////////////////////////

                            Scalar fracture_energy = 0.;
                            for (int d = 0; d < Dim; ++d) {
                                auto gc = grad_x(cell, qp, phase_field * Dim + d);
                                fracture_energy += gc * gc;
                            }

                            fracture_energy *= length_scale / 2.0;
                            auto c = x(cell, qp, phase_field);
                            fracture_energy += 1. / (2.0 * length_scale) * c * c;

                            //////////////////////////////////////////////////////

                            Scalar elastic_energy =
                                (degradation_function.value(cell, qp) * (1.0 - regularization) + regularization) *
                                stress_op.energy(cell, qp);

                            val += elastic_energy + fracture_energy;

                            //////////////////////////////////////////////////////
                            // Pressure induction
                            val +=
                                degradation_function.value(cell, qp) * p(cell, qp) * stress_op.strain_trace(cell, qp);

                            //////////////////////////////////////////////////////
                            integr += val * measure(cell, qp);
                        }

                        data(cell) += integr;
                    });

                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::assemble_scalar");
                return true;
            }

            bool assemble_matrix() override {
                using DegradationFunction = utopia::kokkos::QuadraticDegradationKernel<InterpolateField, Scalar>;

                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_matrix");

                this->ensure_matrix_accumulator();
                auto data = this->matrix_data();

                // FE
                auto &fe = this->fe();
                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int n_var = this->n_vars();
                const int n_quad_points = fe.n_quad_points();

                // Fields
                const int displacement = op_.displacement;
                const int phase_field = op_.phase_field;

                auto x = this->current_solution()->interpolate();
                auto p = pressure_field_->interpolate();
                auto grad_x = gradient_->data();

                // Model parameters
                const Scalar regularization = op_.regularization;
                const Scalar fracture_toughness = op_.fracture_toughness;
                const Scalar length_scale = op_.length_scale;

                // Function and operators
                DegradationFunction degradation_function(x, phase_field);
                CoeffLinearElasticityOp stress_op(op_.lambda, op_.mu, grad_x, measure, displacement);
                LinearElasticityOp elasticity_op(op_.lambda, op_.mu, grad, measure);

                this->loop_cell_test_trial(
                    "IsotropicPhaseFieldForBrittleFractures::assemble_matrix",
                    UTOPIA_LAMBDA(const int cell, const int i, const int j) {
                        const int offset_i = i * n_var;
                        const int offset_j = j * n_var;

                        for (int sub_i = 0; sub_i < Dim; ++sub_i) {
                            int u_dof_i = offset_i + sub_i + displacement;

                            // Integrate displacement block (u, u)
                            for (int sub_j = 0; sub_j < Dim; ++sub_j) {
                                int dof_j = offset_j + sub_j + displacement;

                                Scalar integr = 0.;
                                for (int qp = 0; qp < n_quad_points; ++qp) {
                                    Scalar val = 0.;

                                    val +=
                                        ((1 - regularization) * degradation_function.value(cell, qp) + regularization) *
                                        elasticity_op.inner(cell, i, j, qp, sub_i, sub_j);

                                    integr += val * measure(cell, qp);
                                }

                                data(cell, u_dof_i, dof_j) += integr;
                            }

                            // Integrate coupling block (u, c)
                            int dof_j = offset_j + phase_field;

                            Scalar integr = 0.;
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                Scalar val = 0.;
                                auto ddf = degradation_function.deriv(cell, qp);

                                val += ddf * stress_op.inner_with_strain_test(grad, cell, i, qp, sub_i) * fun(j, qp);

                                //////////////////////////////////////////////////////
                                // Pressure induction
                                val += ddf * p(cell, qp) *
                                       kernels::LinearizedStrain<Dim, Scalar>::trace(grad, cell, i, qp, sub_i) *
                                       fun(j, qp);

                                //////////////////////////////////////////////////////
                                integr += val * measure(cell, qp);
                            }

                            data(cell, u_dof_i, dof_j) += integr;
                        }

                        int u_dof_i = offset_i + phase_field;

                        // Integrate phase_field block (c, c)
                        {
                            int u_dof_j = offset_j + phase_field;

                            Scalar integr = 0.;
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                Scalar val = 0.0;
                                Scalar elastic_energy = stress_op.energy(cell, qp);

                                Scalar lapl = 0.0;
                                for (int d = 0; d < Dim; ++d) {
                                    lapl += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                                }

                                const Scalar diffusion = fracture_toughness * length_scale * lapl;

                                const Scalar reaction = (fracture_toughness / length_scale);

                                const Scalar deriv_cc_elast =
                                    ((1 - regularization) * degradation_function.deriv2(cell, qp)) * elastic_energy;

                                val += diffusion + (reaction + deriv_cc_elast) * fun(i, qp) * fun(j, qp);

                                //////////////////////////////////////////////////////
                                // Pressure induction
                                val += degradation_function.deriv2(cell, qp) * p(cell, qp) *
                                       stress_op.strain_trace(cell, qp) * fun(i, qp) * fun(j, qp);
                                //////////////////////////////////////////////////////

                                integr += val * measure(cell, qp);
                            }
                        }

                        // TODO Integrate coupling block (c, u)

                        for (int sub_j = 0; sub_j < Dim; ++sub_j) {
                            int dof_j = offset_j + sub_j + displacement;

                            Scalar integr = 0.;

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                Scalar val = 0.0;
                                // HERE
                                integr += val * measure(cell, qp);
                            }

                            data(cell, u_dof_i, dof_j) += integr;
                        }
                    });

                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                using DegradationFunction = utopia::kokkos::QuadraticDegradationKernel<InterpolateField, Scalar>;

                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_vector");

                this->ensure_vector_accumulator();
                auto data = this->vector_data();

                // FE
                auto &fe = this->fe();
                auto &&grad = fe.grad();
                auto &&fun = fe.fun();
                auto &&measure = fe.measure();

                const int n_shape_functions = fe.n_shape_functions();
                const int n_var = this->n_vars();
                const int n_quad_points = fe.n_quad_points();

                // Fields
                const int displacement = op_.displacement;
                const int phase_field = op_.phase_field;

                auto x = this->current_solution()->interpolate();
                auto p = pressure_field_->interpolate();
                auto grad_x = gradient_->data();

                // Model parameters
                const Scalar regularization = op_.regularization;
                const Scalar fracture_toughness = op_.fracture_toughness;
                const Scalar length_scale = op_.length_scale;

                // Function and operators
                DegradationFunction degradation_function(x, phase_field);
                CoeffLinearElasticityOp stress_op(op_.lambda, op_.mu, grad_x, measure, displacement);

                this->loop_cell_test(
                    "IsotropicPhaseFieldForBrittleFractures::assemble_vector",
                    UTOPIA_LAMBDA(const int &cell, const int &i) {
                        // Displacement block
                        for (int j = 0; j < n_shape_functions; ++j) {
                            for (int sub_i = 0; sub_i < Dim; ++sub_i) {
                                Scalar integr = 0.;
                                for (int qp = 0; qp < n_quad_points; ++qp) {
                                    Scalar df = degradation_function.value(cell, qp);

                                    Scalar stress_dot_stress =
                                        stress_op.inner_with_strain_test(grad, cell, i, qp, sub_i) *
                                        (df * (1 - regularization) + regularization);

                                    Scalar val = stress_dot_stress;

                                    //////////////////////////////////////////////////////
                                    // Pressure induction
                                    val += df * p(cell, qp) *
                                           kernels::LinearizedStrain<Dim, Scalar>::trace(grad, cell, i, qp, sub_i);
                                    //////////////////////////////////////////////////////

                                    integr += val * measure(cell, qp);
                                }

                                data(cell, i * n_var + displacement + sub_i) += integr;
                            }
                        }

                        // Phase-field block
                        for (int j = 0; j < n_shape_functions; ++j) {
                            Scalar integr = 0.;
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                Scalar val = 0.;
                                Scalar ddf = degradation_function.deriv(cell, qp);
                                Scalar grad_c_elastic_energy =
                                    ddf * (1 - regularization) * stress_op.energy(cell, qp) * fun(i, qp);

                                Scalar grad_c_fracture_energy = 0;
                                for (int d = 0; d < Dim; ++d) {
                                    grad_c_fracture_energy +=
                                        grad_x(cell, qp, phase_field * Dim + d) * grad(cell, i, qp, d);
                                }

                                grad_c_fracture_energy *= length_scale;
                                grad_c_fracture_energy += (1. / length_scale * x(cell, qp, phase_field) * fun(i, qp));
                                grad_c_fracture_energy *= fracture_toughness;

                                val += grad_c_elastic_energy + grad_c_fracture_energy;

                                //////////////////////////////////////////////////////
                                // Pressure induction
                                val += ddf * p(cell, qp) * stress_op.strain_trace(cell, qp) * fun(i, qp);
                                //////////////////////////////////////////////////////

                                integr += val * measure(cell, qp);
                            }

                            data(cell, i * n_var + phase_field) += integr;
                        }
                    });

                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::assemble_vector");
                return true;
            }

            bool apply(const VectorView &x, VectorView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::apply");
                assert(false && "IMPLEMENT ME");
                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::apply");
                return true;
            }

            // // NVCC_PRIVATE :
            Params op_;
            std::shared_ptr<Field<FE>> pressure_field_;
            std::shared_ptr<utopia::kokkos::Gradient<FE>> gradient_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP