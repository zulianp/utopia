#ifndef UTOPIA_INTREPID2_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP
#define UTOPIA_INTREPID2_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP

#include "utopia_Views.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"
#include "utopia_intrepid2_Gradient.hpp"
#include "utopia_intrepid2_PhaseFieldKernels.hpp"
#include "utopia_intrepid2_Strain.hpp"

namespace utopia {

    template <typename Scalar, int Dim>
    class IsotropicPhaseFieldForBrittleFractures : public Configurable {
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
        }

        IsotropicPhaseFieldForBrittleFractures()
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
    };

    namespace intrepid2 {

        template <typename Scalar, int Dim>
        class Assemble<IsotropicPhaseFieldForBrittleFractures<Scalar, Dim>, Scalar> : public FEAssembler<Scalar> {
        public:
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            using UserOp = utopia::IsotropicPhaseFieldForBrittleFractures<Scalar, Dim>;
            using QuadraticDegradation = utopia::intrepid2::QuadraticDegradation<UserOp>;
            using QuadraticPhaseFieldKernels = utopia::intrepid2::PhaseFieldKernels<UserOp, QuadraticDegradation>;

            using V = StaticVector<Scalar, Dim>;
            using M = StaticMatrix<Scalar, Dim, Dim>;
            using SIMDType = Scalar;

            using InterpolateField = typename utopia::intrepid2::Field<Scalar>::Interpolate;
            using InterpolateStrain = typename utopia::intrepid2::LinearizedStrain<Scalar, Dim>::Interpolate;
            using InterpolatePhaseFieldGradient = typename utopia::intrepid2::Gradient<Scalar>::Rank1Op;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim + 1; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return true; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "IsotropicPhaseFieldForBrittleFractures"; }

            virtual bool update(const std::shared_ptr<Field<Scalar>> &displacement) override {
                if (!Super::update(displacement)) {
                    return false;
                }

                //     assert(displacement);
                //     assert(displacement->is_coefficient());

                //     if (!displacement->is_coefficient()) {
                //         Utopia::Abort(
                //             "Assemble<IsotropicPhaseFieldForBrittleFractures>::update, displacement must me in
                //             coefficient " "form!");
                //     }

                //     if (!deformation_gradient_) {
                //         // Initialize gradient
                //         deformation_gradient_ = std::make_shared<Gradient<Scalar>>(this->fe_ptr());
                //     }

                //     deformation_gradient_->init(*displacement);
                //     deformation_gradient_->add_identity();
                //     assert(deformation_gradient_->check_dets_are_positive());
                return true;
            }

            class EnergyIrreversibility {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    Scalar diff = phase_field_(cell, qp) - phase_field_old_(cell, qp);
                    return penalty_param_ / 2 * diff * diff * (diff < 0.0);
                }

                EnergyIrreversibility(const Scalar &penalty_param,
                                      const InterpolateField &phase_field_,
                                      const InterpolateField &phase_field_old_)
                    : penalty_param_(penalty_param), phase_field_(phase_field_), phase_field_old_(phase_field_old_) {}

                Scalar penalty_param_;
                InterpolateField phase_field_;
                InterpolateField phase_field_old_;
            };

            class GradientIrreversibility {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int qp) const {
                    const Scalar diff = phase_field_(cell, qp) - phase_field_old_(cell, qp);
                    return penalty_param_ * diff * (diff < 0.0) * fun_(i, qp);
                }

                GradientIrreversibility(const DynRankView &fun,
                                        const Scalar &penalty_param,
                                        const InterpolateField &phase_field_,
                                        const InterpolateField &phase_field_old_)
                    : fun_(fun),
                      penalty_param_(penalty_param),
                      phase_field_(phase_field_),
                      phase_field_old_(phase_field_old_) {}

                DynRankView fun_;
                Scalar penalty_param_;
                InterpolateField phase_field_;
                InterpolateField phase_field_old_;
            };

            template <class DegradationFunction>
            class Energy {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    const Scalar strain_trace = strain_.trace(cell, qp);
                    const Scalar phase_field_value = phase_field_(cell, qp);

                    Scalar energy = 0.0;
                    if (op_.use_pressure) {
                        energy =
                            DegradationFunction::value(op_, phase_field_value) * pressure_(cell, qp) * strain_trace;
                    }

                    energy += elastic_energy(cell, qp, phase_field_value, strain_trace);
                    energy += fracture_energy(cell, qp, phase_field_value);
                    return energy;
                }

                // UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {

                // }

                UTOPIA_INLINE_FUNCTION Scalar elastic_energy(const int cell,
                                                             const int qp,
                                                             const Scalar &phase_field_value,
                                                             const Scalar &strain_trace) const {
                    return (DegradationFunction::value(op_, phase_field_value) * (1.0 - op_.regularization) +
                            op_.regularization) *
                           strain_energy(strain_trace);
                }

                UTOPIA_INLINE_FUNCTION Scalar strain_energy(const int cell,
                                                            const int qp,
                                                            const Scalar &strain_trace) const {
                    return 0.5 * op_.lambda * strain_trace * strain_trace + op_.mu * strain_.squared_norm(cell, qp);
                }

                UTOPIA_INLINE_FUNCTION Scalar fracture_energy(const int cell,
                                                              const int qp,
                                                              const Scalar &phase_field_value) {
                    return op_.fracture_toughness *
                           (1. / (2.0 * op_.length_scale) * phase_field_value * phase_field_value +
                            op_.length_scale / 2.0 * phase_field_gradient_.squared_norm(cell, qp));
                }

                Energy(const UserOp &op,
                       const InterpolateField &phase_field,
                       const InterpolatePhaseFieldGradient &phase_field_gradient,
                       const InterpolateField &pressure,
                       const InterpolateStrain &strain)
                    : op_(op),
                      phase_field_(phase_field),
                      phase_field_gradient_(phase_field_gradient),
                      pressure_(pressure),
                      strain_(strain) {}

                UserOp op_;
                InterpolateField phase_field_;
                InterpolatePhaseFieldGradient phase_field_gradient_;
                InterpolateField pressure_;
                InterpolateStrain strain_;
            };

            template <class DegradationFunction>
            class ScalarOp {
            public:
                UTOPIA_INLINE_FUNCTION ScalarOp(const Energy<DegradationFunction> &energy, const DynRankView &measure)
                    : energy_(energy), measure_(measure) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    return energy_(cell, qp) * measure_(cell, qp);
                }

                Energy<DegradationFunction> energy_;
                DynRankView measure_;
            };

            template <class DegradationFunction>
            class ScalarOpWithPenalty {
            public:
                UTOPIA_INLINE_FUNCTION ScalarOpWithPenalty(const Energy<DegradationFunction> &energy,
                                                           const DynRankView &measure,
                                                           const EnergyIrreversibility &penalty)
                    : energy_(energy), penalty_(penalty), measure_(measure) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp) const {
                    return energy_(cell, qp) * penalty_(cell, qp) * measure_(cell, qp);
                }

                Energy<DegradationFunction> energy_;
                EnergyIrreversibility penalty_;
                DynRankView measure_;
            };

            class VectorOp {};

            class MatrixOp {};

            // bool assemble_scalar() override {}

            // // FGF_i = F_inv_t_qp * transpose(Grad_i) * F_inv_t_qp
            // UTOPIA_INLINE_FUNCTION static void selective_triple_product(const int di,
            //                                                             const M &F_inv_t_qp,
            //                                                             const V &grad,
            //                                                             M &FGF_i) {
            //     // 2x slower
            //     // M grad_mat;
            //     // make_grad(di, grad, grad_mat);
            //     // FGF_i = F_inv_t_qp * transpose(grad_mat) * F_inv_t_qp;

            //     M GtF;
            //     for (int k1 = 0; k1 < Dim; ++k1) {
            //         for (int k2 = 0; k2 < Dim; ++k2) {
            //             GtF(k1, k2) = grad[k1] * F_inv_t_qp(di, k2);
            //         }
            //     }

            //     FGF_i = F_inv_t_qp * GtF;
            // }

            // UTOPIA_INLINE_FUNCTION static Scalar selective_inner(const int di, const M &mat, const V &grad) {
            //     Scalar ret = 0.0;

            //     for (int dj = 0; dj < Dim; ++dj) {
            //         ret += mat(di, dj) * grad[dj];
            //     }

            //     // 2x slower
            //     // M grad_mat;
            //     // make_grad(di, grad, grad_mat);
            //     // ret = inner(mat, grad_mat);
            //     return ret;
            // }

            // UTOPIA_INLINE_FUNCTION static void make_grad(const int di, const V &vec_grad, M &mat_grad) {
            //     mat_grad.set(0.0);

            //     for (int i = 0; i < Dim; ++i) {
            //         mat_grad(di, i) = vec_grad[i];
            //     }
            // }

            // UTOPIA_INLINE_FUNCTION static void copy_def_grad(const int cell,
            //                                                  const int qp,
            //                                                  const DynRankView &F,
            //                                                  M &F_qp) {
            //     for (int d1 = 0; d1 < Dim; ++d1) {
            //         for (int d2 = 0; d2 < Dim; ++d2) {
            //             F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
            //         }
            //     }
            // }

            // UTOPIA_INLINE_FUNCTION static void linearized_stress(const int sub_trial,
            //                                                      const Scalar &lambda,
            //                                                      const Scalar &mu,
            //                                                      const Scalar &rescale,
            //                                                      const V &grad_trial,
            //                                                      const M &F,
            //                                                      const M &F_inv_t,
            //                                                      const Scalar &J,
            //                                                      M &stress_lin) {
            //     M FGF;
            //     assert(J > 0.0);
            //     assert(J != 0.0);
            //     assert(J == J);

            //     Scalar log_J = device::log(J);
            //     assert(log_J == log_J);

            //     Scalar beta = (rescale * (lambda * log_J - mu));

            //     // FGF = F_inv_t * transpose(Grad_i) * F_inv_t
            //     selective_triple_product(sub_trial, F_inv_t, grad_trial, FGF);
            //     Scalar inner_F_grad_trial = selective_inner(sub_trial, F_inv_t, grad_trial);

            //     make_grad(sub_trial, grad_trial, stress_lin);
            //     stress_lin *= rescale * mu;
            //     stress_lin += (-beta) * FGF + (lambda * inner_F_grad_trial) * F_inv_t;
            // }

            // class Op {
            // public:
            //     Op(const FirstLameParameter &lambda,
            //        const ShearModulus &mu,
            //        const Scalar &rescale,
            //        const DynRankView &F,
            //        const DynRankView &grad,
            //        const DynRankView &measure)
            //         : lambda(lambda),
            //           mu(mu),
            //           rescale(rescale),
            //           F(F),
            //           grad(grad),
            //           measure(measure),
            //           num_qp(measure.extent(1)) {}

            //     UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

            //     UTOPIA_INLINE_FUNCTION Scalar
            //     operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
            //         Scalar ret = 0.0;
            //         for (int qp = 0; qp < num_qp; ++qp) {
            //             ret += (*this)(cell, i, j, qp, sub_i, sub_j);
            //         }

            //         return ret;
            //     }

            //     UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
            //                                              const int i,
            //                                              const int j,
            //                                              const int qp,
            //                                              const int sub_i,
            //                                              const int sub_j) const {
            //         V grad_i, grad_j;
            //         M stress_lin;
            //         M F_qp, F_inv_t_qp;

            //         Scalar dX = measure(cell, qp);

            //         for (int d = 0; d < Dim; ++d) {
            //             grad_i[d] = grad(cell, i, qp, d);
            //             grad_j[d] = grad(cell, j, qp, d);
            //         }

            //         copy_def_grad(cell, qp, F, F_qp);

            //         F_inv_t_qp = inv(transpose(F_qp));
            //         Scalar J = det(F_qp);
            //         linearized_stress(sub_i, lambda, mu, rescale, grad_i, F_qp, F_inv_t_qp, J, stress_lin);
            //         Scalar val = selective_inner(sub_j, stress_lin, grad_j) * dX;
            //         return val;
            //     }

            //     const FirstLameParameter lambda;
            //     const ShearModulus mu;
            //     const Scalar rescale;
            //     const DynRankView F;
            //     const DynRankView grad;
            //     const DynRankView measure;
            //     const int num_qp;
            // };

            // class OpAndStoreHessian {
            // public:
            //     OpAndStoreHessian(const FirstLameParameter &lambda,
            //                       const ShearModulus &mu,
            //                       const Scalar &rescale,
            //                       const DynRankView &F,
            //                       const DynRankView &grad,
            //                       const DynRankView &measure,
            //                       DynRankView &data)
            //         : lambda(lambda),
            //           mu(mu),
            //           rescale(rescale),
            //           F(F),
            //           grad(grad),
            //           measure(measure),
            //           data(data),
            //           num_qp(measure.extent(1)) {}

            //     UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i, const int j) const {
            //         for (int qp = 0; qp < num_qp; ++qp) {
            //             (*this)(cell, i, j, qp);
            //         }
            //     }

            //     UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i, const int j, const int qp) const
            //     {
            //         V grad_i, grad_j;
            //         M FGF_i, stress_lin_i;
            //         M F_qp, F_inv_t_qp;

            //         auto dX = measure(cell, qp);

            //         for (int d = 0; d < Dim; ++d) {
            //             grad_i[d] = grad(cell, i, qp, d);
            //             grad_j[d] = grad(cell, j, qp, d);
            //         }

            //         copy_def_grad(cell, qp, F, F_qp);

            //         F_inv_t_qp = inv(transpose(F_qp));
            //         Scalar J = det(F_qp);
            //         assert(J > 0.0);
            //         assert(J != 0.0);
            //         assert(J == J);

            //         Scalar log_J = device::log(J);
            //         assert(log_J == log_J);

            //         Scalar beta = (rescale * (lambda * log_J - mu));

            //         for (int di = 0; di < Dim; ++di) {
            //             auto dof_i = i * Dim + di;

            //             // FGF_i = F_inv_t_qp * transpose(Grad_i) * F_inv_t_qp
            //             selective_triple_product(di, F_inv_t_qp, grad_i, FGF_i);
            //             Scalar inner_F_grad_i = selective_inner(di, F_inv_t_qp, grad_i);

            //             make_grad(di, grad_i, stress_lin_i);
            //             stress_lin_i *= rescale * mu;
            //             stress_lin_i += (-beta) * FGF_i + (lambda * inner_F_grad_i) * F_inv_t_qp;

            //             for (int dj = 0; dj < Dim; ++dj) {
            //                 auto dof_j = j * Dim + dj;
            //                 auto val = selective_inner(dj, stress_lin_i, grad_j) * dX;
            //                 data(cell, dof_i, dof_j) += val;
            //             }
            //         }
            //     }

            //     const FirstLameParameter lambda;
            //     const ShearModulus mu;
            //     const Scalar rescale;
            //     const DynRankView F;
            //     const DynRankView grad;
            //     const DynRankView measure;
            //     DynRankView data;
            //     const int num_qp;
            // };

            // class OpAndStoreGradient {
            // public:
            //     OpAndStoreGradient(const FirstLameParameter &lambda,
            //                        const ShearModulus &mu,
            //                        const Scalar &rescale,
            //                        const DynRankView &F,
            //                        const DynRankView &grad,
            //                        const DynRankView &measure,
            //                        DynRankView &data)
            //         : lambda(lambda),
            //           mu(mu),
            //           rescale(rescale),
            //           F(F),
            //           grad(grad),
            //           measure(measure),
            //           data(data),
            //           num_qp(measure.extent(1)) {}

            //     UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int j) const {
            //         for (int qp = 0; qp < num_qp; ++qp) {
            //             (*this)(cell, j, qp);
            //         }
            //     }

            //     UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int j, const int qp) const {
            //         V grad_j;
            //         M P;  // P = F_qp
            //         M F_inv_t_qp;

            //         auto dX = measure(cell, qp);

            //         for (int d = 0; d < Dim; ++d) {
            //             grad_j[d] = grad(cell, j, qp, d);
            //         }

            //         copy_def_grad(cell, qp, F, P);

            //         F_inv_t_qp = inv(transpose(P));
            //         Scalar J = det(P);
            //         Scalar log_J = device::log(J);
            //         Scalar lambda_log_J = (rescale * lambda) * log_J;

            //         assert(J > 0.0);
            //         assert(J == J);
            //         assert(log_J == log_J);

            //         P -= F_inv_t_qp;
            //         P *= rescale * mu;
            //         P += lambda_log_J * F_inv_t_qp;

            //         for (int dj = 0; dj < Dim; ++dj) {
            //             auto dof_j = j * Dim + dj;
            //             auto val = selective_inner(dj, P, grad_j) * dX;
            //             data(cell, dof_j) += val;
            //         }
            //     }

            //     const FirstLameParameter lambda;
            //     const ShearModulus mu;
            //     const Scalar rescale;
            //     const DynRankView F;
            //     const DynRankView grad;
            //     const DynRankView measure;
            //     DynRankView data;
            //     const int num_qp;
            // };

            // inline Op make_op() const {
            //     assert(deformation_gradient_);
            //     auto F = deformation_gradient_->data();

            //     return Op(op_.lambda, op_.mu, op_.rescale, F, this->fe().grad, this->fe().measure);
            // }

            // bool apply(const DynRankView &x, DynRankView &y) override {
            //     UTOPIA_TRACE_REGION_BEGIN("Assemble<IsotropicPhaseFieldForBrittleFractures>::apply");

            //     this->apply_vector_operator("Assemble<IsotropicPhaseFieldForBrittleFractures>::apply", x, y,
            //     make_op());

            //     UTOPIA_TRACE_REGION_END("Assemble<IsotropicPhaseFieldForBrittleFractures>::apply");
            //     return false;
            // }

            // bool assemble_matrix() override {
            //     UTOPIA_TRACE_REGION_BEGIN("Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_matrix");

            //     this->ensure_matrix_accumulator();

            //     auto &fe = this->fe();
            //     auto data = this->matrix_data();

            //     {
            //         assert(deformation_gradient_);
            //         auto F = deformation_gradient_->data();
            //         this->loop_cell_test_trial(
            //             "Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_matrix",
            //             OpAndStoreHessian(op_.lambda, op_.mu, op_.rescale, F, fe.grad, fe.measure, data));
            //     }

            //     UTOPIA_TRACE_REGION_END("Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_matrix");
            //     return true;
            // }

            // bool assemble_vector() override {
            //     UTOPIA_TRACE_REGION_BEGIN("Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_vector");

            //     this->ensure_vector_accumulator();

            //     auto &fe = this->fe();
            //     auto data = this->vector_data();

            //     {
            //         auto F = deformation_gradient_->data();
            //         this->loop_cell_test(
            //             "Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_vector",
            //             OpAndStoreGradient(op_.lambda, op_.mu, op_.rescale, F, fe.grad, fe.measure, data));
            //     }

            //     UTOPIA_TRACE_REGION_END("Assemble<IsotropicPhaseFieldForBrittleFractures>::assemble_vector");
            //     return true;
            // }

            // // NVCC_PRIVATE :
            UserOp op_;
            std::shared_ptr<Field<Scalar>> pressure_field_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_ISOTROPIC_PHASE_FIELD_FOR_BRITTLE_FRACTURES_HPP