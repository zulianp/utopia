#ifndef UTOPIA_KOKKOS_COEFF_LINEAR_ELASTICITY_OP_HPP
#define UTOPIA_KOKKOS_COEFF_LINEAR_ELASTICITY_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_kokkos_StrainOp.hpp"

#include "utopia_kokkos_LinearElasticityOp.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim,
                      typename Scalar,
                      class FirstLameParameter,
                      class ShearModulus,
                      class Gradient,
                      class Measure>
            class CoeffLinearElasticityOp {
            public:
                UTOPIA_INLINE_FUNCTION Scalar strain(const int cell,
                                                     const int qp,
                                                     const int sub_i,
                                                     const int sub_j) const {
                    return 0.5 * (grad(cell, qp, (displacement + sub_i) * dim() + sub_j) +
                                  grad(cell, qp, (displacement + sub_j) * dim() + sub_i));
                }

                UTOPIA_INLINE_FUNCTION Scalar strain_trace(const int cell, const int qp) const {
                    Scalar ret = 0.0;

                    for (int d = 0; d < dim(); ++d) {
                        ret += grad(cell, qp, (displacement + d) * dim() + d);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar stress(const int cell,
                                                     const int qp,
                                                     const int sub_i,
                                                     const int sub_j) const {
                    return mux2 * strain(cell, qp, sub_i, sub_j) + lambda * strain_trace(cell, qp);
                }

                UTOPIA_INLINE_FUNCTION Scalar energy(const int cell, const int qp) {
                    Scalar trace = 0.0;
                    Scalar inner_strain = 0.0;
                    for (int d1 = 0; d1 < dim(); ++d1) {
                        trace += strain(cell, qp, d1, d1);

                        for (int d2 = 0; d2 < dim(); ++d2) {
                            const Scalar s2 = strain(cell, qp, d1, d2);
                            inner_strain += s2 * s2;
                        }
                    }

                    return 0.5 * (lambda * trace * trace + mux2 * inner_strain);
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                UTOPIA_INLINE_FUNCTION CoeffLinearElasticityOp(const FirstLameParameter &lambda,
                                                               const ShearModulus &mu,
                                                               const Gradient &grad,
                                                               const Measure &measure,
                                                               const int displacement)
                    : lambda(lambda),
                      mux2(mu * 2),
                      grad(grad),
                      measure(measure),
                      n_qp(measure.extent(1)),
                      displacement(displacement) {}

                template <class TestGrad>
                UTOPIA_INLINE_FUNCTION Scalar inner_with_strain_test(const TestGrad &grad,
                                                                     const int cell,
                                                                     const int i,
                                                                     const int qp,
                                                                     const int sub_i) {
                    Scalar ret = 0.0;

                    for (int d = 0; d < dim(); ++d) {
                        ret += stress(cell, qp, sub_i * dim() + d) * grad(cell, i, qp, d);
                        ret += stress(cell, qp, sub_i + d * dim()) * grad(cell, i, qp, d);
                    }

                    return 0.5 * ret;
                }

                const FirstLameParameter lambda;
                const ShearModulus mux2;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
                // Offset in coefficient array
                const int displacement;
            };
        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia
#endif  // UTOPIA_KOKKOS_COEFF_LINEAR_ELASTICITY_OP_HPP