#ifndef UTOPIA_KOKKOS_GRADIENT_OP_HPP
#define UTOPIA_KOKKOS_GRADIENT_OP_HPP

#include "utopia_Base.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <typename Scalar, class Gradient, class Coefficient>
            class GradientOp {
            public:
                class Rank1 {
                public:
                    UTOPIA_INLINE_FUNCTION Rank1(const Gradient &grad, const Coefficient &coeff)
                        : grad(grad), coeff(coeff), n_shape_functions(grad.extent(1)), n_var(grad.extent(3)) {}

                    UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int d) const {
                        assert(d < n_var);

                        Scalar ret = 0.0;
                        for (int i = 0; i < n_shape_functions; ++i) {
                            ret += coeff(cell, i) * grad(cell, i, qp, d);
                        }

                        return ret;
                    }

                    UTOPIA_INLINE_FUNCTION Scalar squared_norm(const int cell, const int qp) const {
                        Scalar ret = 0.0;

                        for (int d = 0; d < n_var; ++d) {
                            auto x = (*this)(cell, qp, d);
                            ret += x * x;
                        }
                    }

                    UTOPIA_INLINE_FUNCTION int dim() const { return grad.extent(3); }

                    const Gradient grad;
                    const Coefficient coeff;
                    const int n_shape_functions;
                    const int n_var;
                };

                class Rank2 {
                public:
                    UTOPIA_INLINE_FUNCTION Rank2(const Gradient &grad, const Coefficient &coeff)
                        : grad(grad),
                          coeff(coeff),
                          n_shape_functions(grad.extent(1)),
                          n_var(coeff.extent(1) / n_shape_functions) {}

                    UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                             const int qp,
                                                             const int var,
                                                             const int d) const {
                        Scalar ret = 0.0;
                        for (int i = 0; i < n_shape_functions; ++i) {
                            ret += coeff(cell, i * n_var + var) * grad(cell, i, qp, d);
                        }

                        return ret;
                    }

                    UTOPIA_INLINE_FUNCTION int dim() const { return grad.extent(3); }

                    const Gradient grad;
                    const Coefficient coeff;
                    const int n_shape_functions;
                    const int n_var;
                };
            };

            template <class Op, class OuputField>
            class StoreRank1Gradient {
            public:
                UTOPIA_INLINE_FUNCTION StoreRank1Gradient(const Op &op, OuputField &field) : op(op), field(field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int d) const {
                    field(cell, qp, d) = op(cell, qp, d);
                }

                typename Op::Rank1 op;
                OuputField field;
            };

            template <class Op, class OuputField>
            class StoreRank2Gradient {
            public:
                UTOPIA_INLINE_FUNCTION StoreRank2Gradient(const Op &op, OuputField &field)
                    : op(op), spatial_dim(op.grad.extent(3)), field(field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp, const int var, const int d) const {
                    field(cell, qp, var * spatial_dim + d) = op(cell, qp, var, d);
                }

                typename Op::Rank2 op;
                const int spatial_dim;
                OuputField field;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_GRADIENT_OP_HPP
