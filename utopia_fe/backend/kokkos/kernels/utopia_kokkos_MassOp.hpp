#ifndef UTOPIA_KOKKOS_MASS_OP_HPP
#define UTOPIA_KOKKOS_MASS_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <typename Scalar_, class Density, class Fun, class Measure>
            class MassOp {
            public:
                // To be used with LumedOp
                using Scalar = Scalar_;

                UTOPIA_INLINE_FUNCTION MassOp(const Density &density,
                                              const Fun &fun,
                                              const Measure &measure,
                                              const int n_components)
                    : density(density),
                      fun(fun),
                      measure(measure),
                      n_components(n_components),
                      n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j) const {
                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);
                        ret += fun(i, qp) * fun(j, qp) * density * dX;
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j, const int qp) const {
                    auto dX = measure(cell, qp);
                    return fun(i, qp) * fun(j, qp) * density * dX;
                }

                /////////////////////////////////////////////////////////////////////////////////////////////////

                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);
                        ret += fun(i, qp) * fun(j, qp) * density * dX;
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int i,
                                                         const int j,
                                                         const int sub_i,
                                                         const int sub_j,
                                                         const int qp) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    auto dX = measure(cell, qp);
                    return fun(i, qp) * fun(j, qp) * density * dX;
                }

                inline int dim() const { return n_components; }

                const Density density;
                const Fun fun;
                const Measure measure;
                const int n_components;
                const int n_qp;
            };

            template <class WrappedOp>
            class LumpedOp {
            public:
                using Scalar = typename WrappedOp::Scalar;

                UTOPIA_INLINE_FUNCTION LumpedOp(WrappedOp op) : op_(op), n_shape_functions(op.fun.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i) const {
                    Scalar ret = 0.0;
                    for (int j = 0; j < n_shape_functions; ++j) {
                        ret += op_(cell, i, j);
                    }
                    return ret;
                }

                int dim() const { return op_.dim(); }

                WrappedOp op_;
                const int n_shape_functions;
            };

            template <class ElementMatrices>
            class Lump {
            public:
                using Scalar = typename Traits<ElementMatrices>::Scalar;

                UTOPIA_INLINE_FUNCTION Lump(const ElementMatrices &element_matrices)
                    : element_matrices(element_matrices), n_shape_functions(element_matrices.extent(2)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j) const {
                    if (i == j) {
                        return compute(cell, i);
                    } else {
                        return 0.0;
                    }
                }

                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    if (i == j) {
                        return compute(cell, i);
                    } else {
                        return 0.0;
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i) const {
                    auto d = compute(cell, i);

                    for (int j = 0; j < n_shape_functions; ++j) {
                        element_matrices(cell, i, j) = 0.0;
                    }

                    element_matrices(cell, i, i) = d;
                }

                UTOPIA_INLINE_FUNCTION Scalar compute(const int cell, const int i) const {
                    Scalar ret = 0.0;

                    for (int j = 0; j < n_shape_functions; ++j) {
                        ret += element_matrices(cell, i, j);
                    }

                    return ret;
                }

                ElementMatrices element_matrices;
                const int n_shape_functions;
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MASS_OP_HPP
