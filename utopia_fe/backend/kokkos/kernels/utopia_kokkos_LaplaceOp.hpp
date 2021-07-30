#ifndef UTOPIA_KOKKOS_LAPLACE_OP_HPP
#define UTOPIA_KOKKOS_LAPLACE_OP_HPP

#include "utopia_Base.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <typename Scalar, class Coeff, class Gradient, class Measure>
            class LaplaceOp {
            public:
                UTOPIA_INLINE_FUNCTION LaplaceOp(const Coeff &coeff, const Gradient &grad, const Measure &measure)
                    : coeff(coeff), grad(grad), measure(measure), n_qp(measure.extent(1)), dim(grad.extent(3)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar dot_g = 0.0;
                        for (int d = 0; d < dim; ++d) {
                            dot_g += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += dot_g * measure(cell, qp);
                    }

                    return ret * coeff;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell,
                                                         const int &i,
                                                         const int &j,
                                                         const int qp) const {
                    Scalar ret = 0.0;

                    for (int d = 0; d < dim; ++d) {
                        ret += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                    }

                    return ret * coeff * measure(cell, qp);
                }

                const Coeff coeff;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
                const int dim;
            };

            template <int Dim, typename Scalar, class Coeff, class Gradient, class Measure>
            class VectorLaplaceOp {
            public:
                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;

                        for (int d = 0; d < dim(); ++d) {
                            val += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += coeff * val * measure(cell, qp);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                UTOPIA_INLINE_FUNCTION VectorLaplaceOp(const Coeff &coeff, const Gradient &grad, const Measure &measure)
                    : coeff(coeff), grad(grad), measure(measure), n_qp(measure.extent(1)) {}

                const Coeff coeff;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LAPLACE_OP_HPP
