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

                GradientIrreversibility(const Function &fun,
                                        const Scalar &penalty_param,
                                        const InterpolateField &phase_field_,
                                        const InterpolateField &phase_field_old_)
                    : fun_(fun),
                      penalty_param_(penalty_param),
                      phase_field_(phase_field_),
                      phase_field_old_(phase_field_old_) {}

                Function fun_;
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

                Energy(const Params &op,
                       const InterpolateField &phase_field,
                       const InterpolatePhaseFieldGradient &phase_field_gradient,
                       const InterpolateField &pressure,
                       const InterpolateStrain &strain)
                    : op_(op),
                      phase_field_(phase_field),
                      phase_field_gradient_(phase_field_gradient),
                      pressure_(pressure),
                      strain_(strain) {}

                Params op_;
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

            bool assemble_scalar() override {
                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_scalar");

                this->ensure_scalar_accumulator();

                auto &fe = this->fe();
                auto data = this->scalar_data();

                assert(false && "IMPLEMENT ME");
                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::assemble_scalar");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                assert(false && "IMPLEMENT ME");
                UTOPIA_TRACE_REGION_END("IsotropicPhaseFieldForBrittleFractures::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("IsotropicPhaseFieldForBrittleFractures::assemble_vector");

                this->ensure_vector_accumulator();

                auto &fe = this->fe();
                auto data = this->vector_data();
                assert(false && "IMPLEMENT ME");
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