#ifndef UTOPIA_UNIFORM_FE_HPP
#define UTOPIA_UNIFORM_FE_HPP

#include "utopia_Base.hpp"
#include "utopia_kokkos_FEBase.hpp"

#include <Kokkos_DynRankView.hpp>

#include <cassert>

namespace utopia {

    namespace kokkos {

        // template <typename T>
        // using DynView = ::Kokkos::DynRankView<T, ::Kokkos::DefaultExecutionSpace>;

        // template <typename T>
        // using View = ::Kokkos::View<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename Scalar, class ExecutionSpace_ = ::Kokkos::DefaultExecutionSpace>
        class UniformFEDefaultTraits {
        public:
            using ExecutionSpace = ExecutionSpace_;
            using SizeType = int;
            using MeasureView = ::Kokkos::View<Scalar *, ExecutionSpace_>;
            using FunctionView = ::Kokkos::View<Scalar **, ExecutionSpace_>;
            using GradientView = ::Kokkos::View<Scalar ***, ExecutionSpace_>;
            using JacobianView = ::Kokkos::View<Scalar ***, ExecutionSpace_>;
            using JacobianInverseView = ::Kokkos::View<Scalar ***, ExecutionSpace_>;
            using IntView = ::Kokkos::View<int *, ExecutionSpace_>;
            using DynRankView = ::Kokkos::DynRankView<Scalar, ExecutionSpace_>;
        };

        template <typename Scalar_, class Traits = UniformFEDefaultTraits<Scalar_>>
        class UniformFE : public FEBase<typename Traits::ExecutionSpace, typename Traits::SizeType> {
        public:
            using Scalar = Scalar_;
            using SizeType = typename Traits::SizeType;
            using MeasureView = typename Traits::MeasureView;
            using FunctionView = typename Traits::FunctionView;
            using GradientView = typename Traits::GradientView;
            using JacobianView = typename Traits::JacobianView;
            using JacobianInverseView = typename Traits::JacobianInverseView;
            using IntView = typename Traits::IntView;
            using DynRankView = typename Traits::DynRankView;

            virtual ~UniformFE() = default;

            class Gradient {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType fun,
                                                         SizeType qp,
                                                         SizeType d) const {
                    assert(fun < SizeType(grad.extent(0)));
                    assert(qp < SizeType(grad.extent(1)));
                    assert(d < SizeType(grad.extent(2)));

                    return grad(fun, qp, d);
                }

                UTOPIA_INLINE_FUNCTION SizeType extent(const SizeType i) const {
                    if (i == 0) {
                        return n_cells;
                    } else {
                        return grad.extent(i - 1);
                    }
                }

                GradientView grad;
                SizeType n_cells;
            };

            class Function {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType i, SizeType qp) const {
                    assert(i < SizeType(fun.extent(0)));
                    assert(qp < SizeType(fun.extent(1)));

                    return fun(i, qp);
                }

                UTOPIA_INLINE_FUNCTION SizeType extent(const SizeType i) const { return fun.extent(i); }

                FunctionView fun;
            };

            class Jacobian {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType row,
                                                         SizeType col,
                                                         SizeType qp) const {
                    assert(row < SizeType(jacobian.extent(0)));
                    assert(col < SizeType(jacobian.extent(1)));
                    assert(qp < SizeType(jacobian.extent(2)));
                    return jacobian(row, col, qp);
                }

                UTOPIA_INLINE_FUNCTION SizeType extent(const SizeType i) const {
                    if (i == 0) {
                        return n_cells;
                    } else {
                        return jacobian.extent(i - 1);
                    }
                }

                JacobianView jacobian;
                SizeType n_cells;
            };

            class JacobianInverse {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/,
                                                         SizeType row,
                                                         SizeType col,
                                                         SizeType qp) const {
                    assert(row < SizeType(jacobian_inverse.extent(0)));
                    assert(col < SizeType(jacobian_inverse.extent(1)));
                    assert(qp < SizeType(jacobian_inverse.extent(2)));

                    return jacobian_inverse(row, col, qp);
                }

                UTOPIA_INLINE_FUNCTION SizeType extent(const SizeType i) const {
                    if (i == 0) {
                        return n_cells;
                    } else {
                        return jacobian_inverse.extent(i - 1);
                    }
                }

                JacobianInverseView jacobian_inverse;
                SizeType n_cells;
            };

            class Measure {
            public:
                UTOPIA_INLINE_FUNCTION Scalar operator()(SizeType /*cell*/, SizeType qp) const {
                    assert(qp < SizeType(measure.extent(0)));
                    return measure(qp);
                }

                UTOPIA_INLINE_FUNCTION SizeType extent(const SizeType i) const {
                    if (i == 0) {
                        return n_cells;
                    } else {
                        return measure.extent(i - 1);
                    }
                }

                MeasureView measure;
                SizeType n_cells;
            };

            UTOPIA_INLINE_FUNCTION Gradient grad() const { return {grad_, n_cells_}; }
            UTOPIA_INLINE_FUNCTION Function fun() const { return {fun_}; }
            UTOPIA_INLINE_FUNCTION Measure measure() const { return {measure_, n_cells_}; }
            UTOPIA_INLINE_FUNCTION Jacobian jacobian() const { return {jacobian_, n_cells_}; }
            UTOPIA_INLINE_FUNCTION JacobianInverse jacobian_inverse() const { return {jacobian_inverse_, n_cells_}; }

            UTOPIA_INLINE_FUNCTION IntView &element_tags() { return element_tags_; }
            UTOPIA_INLINE_FUNCTION const IntView &element_tags() const { return element_tags_; }

            void init(const SizeType n_cells,
                      const MeasureView &measure,
                      const FunctionView &fun,
                      const GradientView &grad,
                      const JacobianView &jacobian,
                      const JacobianInverseView &jacobian_inverse) {
                n_cells_ = n_cells;
                measure_ = measure;
                fun_ = fun;
                grad_ = grad;
                jacobian_ = jacobian;
                jacobian_inverse_ = jacobian_inverse;
            }

            inline SizeType n_cells() const override { return n_cells_; }
            inline SizeType n_shape_functions() const override { return fun_.extent(0); }
            inline SizeType n_quad_points() const override { return measure_.extent(0); }
            inline SizeType spatial_dimension() const override { return grad_.extent(2); }
            inline SizeType manifold_dimension() const override { return jacobian_.extent(0); }

        private:
            SizeType n_cells_{0};
            MeasureView measure_;
            FunctionView fun_;
            GradientView grad_;
            JacobianView jacobian_;
            JacobianInverseView jacobian_inverse_;

            IntView element_tags_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_UNIFORM_FE_HPP
