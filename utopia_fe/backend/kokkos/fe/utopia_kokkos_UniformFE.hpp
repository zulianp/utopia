#ifndef UTOPIA_UNIFORM_FE_HPP
#define UTOPIA_UNIFORM_FE_HPP

#include "utopia_Base.hpp"

#include <Kokkos_DynRankView.hpp>
#include <Kokkos_View.hpp>

#include <cassert>

namespace utopia {

    namespace kokkos {

        template <typename T>
        using DynView = ::Kokkos::DynRankView<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename T>
        using View = ::Kokkos::View<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename Scalar, class ExecutionSpace = ::Kokkos::DefaultExecutionSpace>
        class UniformFEDefaultTraits {
        public:
            using SizeType = int;
            using MeasureView = utopia::kokkos::View<Scalar *>;
            using FunctionView = utopia::kokkos::View<Scalar **>;
            using GradientView = utopia::kokkos::View<Scalar ***>;
            using JacobianView = utopia::kokkos::View<Scalar ***>;
            using JacobianInverseView = utopia::kokkos::View<Scalar ***>;
        };

        template <typename Scalar, class Traits = UniformFEDefaultTraits<Scalar>>
        class UniformFE {
        public:
            using SizeType = typename Traits::SizeType;
            using MeasureView = typename Traits::MeasureView;
            using FunctionView = typename Traits::FunctionView;
            using GradientView = typename Traits::GradientView;
            using JacobianView = typename Traits::JacobianView;
            using JacobianInverseView = typename Traits::JacobianInverseView;

            class Gradient {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType fun,
                                                         SizeType qp,
                                                         SizeType d) const {
                    assert(fun < grad.extent(0));
                    assert(qp < grad.extent(1));
                    assert(d < grad.extent(2));

                    return grad(fun, qp, d);
                }

                GradientView grad;
            };

            class Function {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType i, SizeType qp) const {
                    assert(i < fun.extent(0));
                    assert(qp < fun.extent(1));

                    return fun(i, qp);
                }

                FunctionView fun;
            };

            class Jacobian {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType row,
                                                         SizeType col,
                                                         SizeType qp) const {
                    assert(row < jacobian.extent(0));
                    assert(col < jacobian.extent(1));
                    assert(qp < jacobian.extent(2));
                    return jacobian(row, col, qp);
                }

                JacobianView jacobian;
            };

            class JacobianInverse {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType row,
                                                         SizeType col,
                                                         SizeType qp) const {
                    assert(row < jacobian_inverse.extent(0));
                    assert(col < jacobian_inverse.extent(1));
                    assert(qp < jacobian_inverse.extent(2));

                    return jacobian_inverse(row, col, qp);
                }

                JacobianInverseView jacobian_inverse;
            };

            class Measure {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/, SizeType qp) const {
                    assert(qp < measure.extent(0));
                    return measure(qp);
                }

                MeasureView measure;
            };

            UTOPIA_INLINE_FUNCTION Gradient grad() const { return {grad_}; }
            UTOPIA_INLINE_FUNCTION Function fun() const { return {fun_}; }
            UTOPIA_INLINE_FUNCTION Measure measure() const { return {measure_}; }
            UTOPIA_INLINE_FUNCTION Jacobian jacobian() const { return {jacobian_}; }
            UTOPIA_INLINE_FUNCTION JacobianInverse jacobian_inverse() const { return {jacobian_inverse_}; }

            void init(const MeasureView &measure,
                      const FunctionView &fun,
                      const GradientView &grad,
                      const JacobianView &jacobian,
                      const JacobianInverseView &jacobian_inverse) {
                measure_ = measure;
                fun_ = fun;
                grad_ = grad;
                jacobian_ = jacobian;
                jacobian_inverse_ = jacobian_inverse;
            }

            // inline SizeType num_cells() const { return cell_nodes.extent(0); }
            UTOPIA_INLINE_FUNCTION SizeType n_shape_functions() const { return fun_.extent(0); }
            UTOPIA_INLINE_FUNCTION SizeType n_quad_points() const { return measure_.extent(0); }
            UTOPIA_INLINE_FUNCTION SizeType spatial_dimension() const { return grad_.extent(2); }
            // inline SizeType manifold_dimension() const { return type.getDimension(); }

        private:
            // View<Scalar**> q_points_;
            // View<Scalar*> q_weights_;

            MeasureView measure_;
            FunctionView fun_;
            GradientView grad_;
            JacobianView jacobian_;
            JacobianInverseView jacobian_inverse_;

            // View<int*> element_tags_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_UNIFORM_FE_HPP
