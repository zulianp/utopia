#ifndef UTOPIA_FE_HPP
#define UTOPIA_FE_HPP

#include "utopia_Base.hpp"

#include <Kokkos_DynRankView.hpp>
#include <Kokkos_View.hpp>

namespace utopia {

    namespace kokkos {

        template <typename T>
        using DynView = ::Kokkos::DynRankView<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename T>
        using View = ::Kokkos::View<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename Scalar, class ExecutionSpace = ::Kokkos::DefaultExecutionSpace>
        class FEDefaultTraits {
        public:
            using SizeType = int;
            using MeasureView = utopia::kokkos::View<Scalar **>;
            using FunctionView = utopia::kokkos::View<Scalar **>;
            using GradientView = utopia::kokkos::View<Scalar ****>;
            using JacobianView = utopia::kokkos::View<Scalar ****>;
            using JacobianInverseView = utopia::kokkos::View<Scalar ****>;
        };

        template <typename Scalar, class Traits = FEDefaultTraits<Scalar>>
        class FE {
        public:
            using SizeType = typename Traits::SizeType;
            using MeasureView = typename Traits::MeasureView;
            using FunctionView = typename Traits::FunctionView;
            using GradientView = typename Traits::GradientView;
            using JacobianView = typename Traits::JacobianView;
            using JacobianInverseView = typename Traits::JacobianInverseView;

            using Measure = typename Traits::MeasureView;
            using Function = typename Traits::FunctionView;
            using Gradient = typename Traits::GradientView;
            using Jacobian = typename Traits::JacobianView;
            using JacobianInverse = typename Traits::JacobianInverseView;

            UTOPIA_INLINE_FUNCTION GradientView grad() const { return grad_; }
            UTOPIA_INLINE_FUNCTION FunctionView fun() const { return fun_; }
            UTOPIA_INLINE_FUNCTION MeasureView measure() const { return measure_; }
            UTOPIA_INLINE_FUNCTION JacobianView jacobian() const { return jacobian_; }
            UTOPIA_INLINE_FUNCTION JacobianInverseView jacobian_inverse() const { return jacobian_inverse_; }

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

#endif  // UTOPIA_FE_HPP
