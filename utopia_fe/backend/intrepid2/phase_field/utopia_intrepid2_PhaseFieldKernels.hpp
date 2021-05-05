#ifndef UTOPIA_INTREPID2_PHASE_FIELD_KERNELS_HPP
#define UTOPIA_INTREPID2_PHASE_FIELD_KERNELS_HPP

#include "utopia_Views.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class Parameters>
        class QuadraticDegradation {
        public:
            template <typename C>
            UTOPIA_INLINE_FUNCTION static C fun(const Parameters &, const C &c) {
                C imc = 1.0 - c;
                return imc * imc;
            }

            template <typename C>
            UTOPIA_INLINE_FUNCTION static C deriv(const Parameters &, const C &c) {
                C imc = 1.0 - c;
                return -2.0 * imc;
            }

            template <typename C>
            UTOPIA_INLINE_FUNCTION static C deriv2(const Parameters &, const C &) {
                return 2.0;
            }
        };

        template <class Parameters, class DegradationFunction>
        class PhaseFieldKernels {
        public:
            template <typename PhaseFieldValue, class StressShape, class Grad>
            UTOPIA_INLINE_FUNCTION static PhaseFieldValue bilinear_uu(const Parameters &params,
                                                                      const PhaseFieldValue &phase_field_value,
                                                                      const StressShape &stress,
                                                                      const Grad &strain_test) {
                const auto gc = ((1.0 - params.regularization) * DegradationFunction::fun(params, phase_field_value) +
                                 params.regularization);

                return inner(gc * stress, strain_test);
            }

            template <class Grad>
            UTOPIA_INLINE_FUNCTION static auto diffusion_c(const Parameters &params,
                                                           const Grad &g_trial,
                                                           const Grad &g_test) {
                return params.fracture_toughness * params.length_scale * inner(g_trial, g_test);
            }

            template <typename FunValue>
            UTOPIA_INLINE_FUNCTION static FunValue reaction_c(const Parameters &params,
                                                              const FunValue &trial,
                                                              const FunValue &test) {
                return (params.fracture_toughness / params.length_scale) * trial * test;
            }

            template <typename PhaseFieldValue, typename ElasticEnergy, typename TrialFunction, typename TestFunction>
            UTOPIA_INLINE_FUNCTION static auto elastic_deriv_cc(const Parameters &params,
                                                                const PhaseFieldValue &phase_field_value,
                                                                const ElasticEnergy &elastic_energy_positive,
                                                                const TrialFunction &trial,
                                                                const TestFunction &test) {
                const auto dcc = (1.0 - params.regularization) * DegradationFunction::deriv2(params, phase_field_value);
                return dcc * trial * elastic_energy_positive * test;
            }

            template <typename PhaseFieldValue, class Grad, typename TestFunction, class GradTest>
            UTOPIA_INLINE_FUNCTION static PhaseFieldValue grad_fracture_energy_wrt_c(
                const Parameters &params,
                const PhaseFieldValue &phase_field_value,
                const Grad &phase_field_grad,
                const TestFunction &test_function,
                const GradTest &grad_test_function) {
                return params.fracture_toughness *
                       ((1. / params.length_scale * phase_field_value * test_function) +
                        (params.length_scale * inner(phase_field_grad, grad_test_function)));
            }

            template <typename PhaseFieldValue, class Stress, class FullStrain, typename CTrialFunction>
            UTOPIA_INLINE_FUNCTION static auto bilinear_uc(const Parameters &params,
                                                           const PhaseFieldValue &phase_field_value,
                                                           const Stress &stress_p,
                                                           const FullStrain &full_strain,
                                                           const CTrialFunction &c_trial_fun) {
                return c_trial_fun *
                       inner(DegradationFunction::deriv(params, phase_field_value) * stress_p, full_strain);
            }

            template <typename PhaseFieldValue, class Grad>
            UTOPIA_INLINE_FUNCTION static PhaseFieldValue fracture_energy(const Parameters &params,
                                                                          const PhaseFieldValue &phase_field_value,
                                                                          const Grad &phase_field_grad) {
                return params.fracture_toughness *
                       (1. / (2.0 * params.length_scale) * phase_field_value * phase_field_value +
                        params.length_scale / 2.0 * inner(phase_field_grad, phase_field_grad));
            }

            template <typename PhaseFieldValue,
                      typename ElasticEnergy,
                      typename TrialFunction,
                      typename TestFunction,
                      typename TrialGrad,
                      typename TestGrad>
            UTOPIA_INLINE_FUNCTION static auto bilinear_cc(const Parameters &params,
                                                           const PhaseFieldValue &phase_field_value,
                                                           const ElasticEnergy &elastic_energy_p,
                                                           const TrialFunction &shape_trial,
                                                           const TestFunction &shape_test,
                                                           const TrialGrad &grad_trial,
                                                           const TestGrad &grad_test) {
                return diffusion_c(params, grad_trial, grad_test) + reaction_c(params, shape_trial, shape_test) +
                       elastic_deriv_cc(params, phase_field_value, elastic_energy_p, shape_trial, shape_test);
            }
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_PHASE_FIELD_KERNELS_HPP