#ifndef UTOPIA_KOKKOS_LINEAR_ELASTICITY_OP_HPP
#define UTOPIA_KOKKOS_LINEAR_ELASTICITY_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_kokkos_StrainOp.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim,
                      typename Scalar,
                      class FirstLameParameter,
                      class ShearModulus,
                      class Gradient,
                      class Measure>
            class LinearElasticityOp {
            public:
                using StrainKernel = kernels::LinearizedStrain<Dim, Scalar>;

                UTOPIA_INLINE_FUNCTION LinearElasticityOp(const FirstLameParameter &lambda,
                                                          const ShearModulus &mu,
                                                          const Gradient &grad,
                                                          const Measure &measure)
                    : lambda(lambda), mux2(mu * 2), grad(grad), measure(measure), n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar strain_inner(const int cell,
                                                           const int i,
                                                           const int j,
                                                           const int qp,
                                                           const int sub_i,
                                                           const int sub_j) const {
                    auto ret = StrainKernel::inner(grad, cell, i, j, qp, sub_i, sub_j);
                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar strain_trace(const int cell,
                                                           const int i,
                                                           const int qp,
                                                           const int sub_i) const {
                    return StrainKernel::trace(grad, cell, i, qp, sub_i);
                }

                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    Scalar ret = 0.0;

                    for (int qp = 0; qp < n_qp; ++qp) {
                        const Scalar val = stress(cell, i, j, qp, sub_i, sub_j);

                        ret += val * measure(cell, qp);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar
                stress(const int cell, const int i, const int j, const int qp, const int sub_i, const int sub_j) const {
                    return mux2 * strain_inner(cell, i, j, qp, sub_i, sub_j) +
                           lambda * strain_trace(cell, i, qp, sub_i) * strain_trace(cell, j, qp, sub_j);
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                const FirstLameParameter lambda;
                const ShearModulus mux2;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LINEAR_ELASTICITY_OP_HPP
