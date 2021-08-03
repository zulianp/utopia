#ifndef UTOPIA_KOKKOS_STRAIN_OP_HPP
#define UTOPIA_KOKKOS_STRAIN_OP_HPP

#include "utopia_Traits.hpp"

#include "utopia_kokkos_GradientOp.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim, typename Scalar>
            class LinearizedStrain {
            public:
                template <class Grad>
                UTOPIA_INLINE_FUNCTION static constexpr auto inner(const Grad &grad,
                                                                   const int cell,
                                                                   const int i,
                                                                   const int j,
                                                                   const int qp,
                                                                   const int sub_i,
                                                                   const int sub_j) {
                    if (sub_i == sub_j) {
                        Scalar ret = 0.0;
                        for (int d = 0; d < Dim; ++d) {
                            ret += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret *= 0.5;
                        ret += 0.5 * grad(cell, i, qp, sub_i) * grad(cell, j, qp, sub_i);
                        return ret;

                    } else {
                        Scalar ret = 0.5 * grad(cell, i, qp, sub_j) * grad(cell, j, qp, sub_i);
                        return ret;
                    }
                }

                template <class Grad>
                UTOPIA_INLINE_FUNCTION static constexpr auto trace(const Grad &grad,
                                                                   const int cell,
                                                                   const int i,
                                                                   const int qp,
                                                                   const int sub_i) {
                    return grad(cell, i, qp, sub_i);
                }
            };

            template <int Dim, typename Scalar, class Gradient, class Coefficient>
            class InterpolateLinearizedStrainOp {
            public:
                using GradOp = typename GradientOp<Scalar, Gradient, Coefficient>::Rank2;

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int qp,
                                                         const int sub_i,
                                                         const int sub_j) const {
                    return 0.5 * (grad(cell, qp, sub_i, sub_j) + grad(cell, qp, sub_j, sub_i));
                }

                UTOPIA_INLINE_FUNCTION Scalar trace(const int cell, const int qp) const {
                    Scalar ret = (*this)(cell, qp, 0, 0);

                    for (int d = 1; d < Dim; ++d) {
                        ret += (*this)(cell, qp, d, d);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar squared_norm(const int cell, const int qp) const {
                    Scalar ret = 0.0;

                    for (int d1 = 0; d1 < Dim; ++d1) {
                        for (int d2 = 0; d2 < Dim; ++d2) {
                            auto x = (*this)(cell, qp, d1, d2);
                            ret += x * x;
                        }
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION InterpolateLinearizedStrainOp(const Gradient &grad, const Coefficient &coeff)
                    : grad(grad, coeff) {}

                GradOp grad;
            };

            template <class StrainOp, class OutputField>
            class StoreInterpolatedStrain {
            public:
                UTOPIA_INLINE_FUNCTION StoreInterpolatedStrain(const StrainOp &op, OutputField &data)
                    : op(op), data(data), n(op.grad.dim()) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int qp,
                                                       const int sub_i,
                                                       const int sub_j) const {
                    data(cell, qp, sub_i * n + sub_j) = op(cell, qp, sub_i, sub_j);
                }

                StrainOp op;
                OutputField data;
                const int n;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_STRAIN_OP_HPP
