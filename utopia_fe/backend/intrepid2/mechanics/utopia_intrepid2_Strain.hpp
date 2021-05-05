#ifndef UTOPIA_INTREPID2_STRAIN_HPP
#define UTOPIA_INTREPID2_STRAIN_HPP

#include "utopia_Traits.hpp"

#include "utopia_intrepid2_Gradient.hpp"

namespace utopia {
    namespace intrepid2 {

        template <int Dim>
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
                using Scalar = typename Traits<Grad>::Scalar;

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

            template <typename Scalar>
            class Interpolate {
            public:
                using GradOp = typename Gradient<Scalar>::Rank2Op;
                using FE = utopia::intrepid2::FE<Scalar>;
                using SizeType = typename FE::SizeType;
                using DynRankView = typename FE::DynRankView;

                UTOPIA_INLINE_FUNCTION Scalar apply(const int cell,
                                                    const int qp,
                                                    const int sub_i,
                                                    const int sub_j) const {
                    return 0.5 * (grad_(cell, qp, sub_i, sub_j) + grad_(cell, qp, sub_j, sub_i));
                }

                UTOPIA_INLINE_FUNCTION Interpolate(const DynRankView &grad, const DynRankView &coeff)
                    : grad_(grad, coeff) {}

                GradOp grad_;
            };

            // UTOPIA_INLINE_FUNCTION static void make(const int dim,
            //                                         Scalar *grad,
            //                                         StaticMatrix<Scalar, Dim, Dim> &strain) {
            //     strain.set(0.0);

            //     for (int i = 0; i < Dim; ++i) {
            //         strain(dim, i) = grad[i];
            //     }

            //     strain.symmetrize();
            // }
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_STRAIN_HPP