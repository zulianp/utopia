#ifndef UTOPIA_FE_HPP
#define UTOPIA_FE_HPP

#include "utopia_Logger.hpp"

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_kokkos_FEBase.hpp"
#include "utopia_kokkos_FEForwardDeclarations.hpp"

#include <Kokkos_DynRankView.hpp>
#include <Kokkos_View.hpp>

namespace utopia {

    namespace kokkos {

        // template <typename T>
        // using DynView = ::Kokkos::DynRankView<T, ::Kokkos::DefaultExecutionSpace>;

        // template <typename T>
        // using View = ::Kokkos::View<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename Scalar, class ExecutionSpace_ = ::Kokkos::DefaultExecutionSpace>
        class FEDefaultTraits {
        public:
            using ExecutionSpace = ExecutionSpace_;
            using SizeType = int;
            using MeasureView = ::Kokkos::View<Scalar **, ExecutionSpace_>;
            using FunctionView = ::Kokkos::View<Scalar **, ExecutionSpace_>;
            using GradientView = ::Kokkos::View<Scalar ****, ExecutionSpace_>;
            using JacobianView = ::Kokkos::View<Scalar ****, ExecutionSpace_>;
            using JacobianInverseView = ::Kokkos::View<Scalar ****, ExecutionSpace_>;
            using IntView = ::Kokkos::View<int *>;
            using DynRankView = ::Kokkos::DynRankView<Scalar, ExecutionSpace_>;
        };

        template <typename Scalar_, class Traits = FEDefaultTraits<Scalar_>>
        class FE : public FEBase<typename Traits::ExecutionSpace, typename Traits::SizeType> {
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

            using Measure = typename Traits::MeasureView;
            using Function = typename Traits::FunctionView;
            using Gradient = typename Traits::GradientView;
            using Jacobian = typename Traits::JacobianView;
            using JacobianInverse = typename Traits::JacobianInverseView;

            virtual ~FE() = default;

            UTOPIA_INLINE_FUNCTION GradientView grad() const { return grad_; }
            UTOPIA_INLINE_FUNCTION FunctionView fun() const { return fun_; }
            UTOPIA_INLINE_FUNCTION MeasureView measure() const { return measure_; }
            UTOPIA_INLINE_FUNCTION JacobianView jacobian() const { return jacobian_; }
            UTOPIA_INLINE_FUNCTION JacobianInverseView jacobian_inverse() const { return jacobian_inverse_; }

            inline bool empty() const { return n_cells() == 0; }
            inline SizeType n_cells() const override { return measure_.extent(0); }
            inline SizeType n_shape_functions() const override { return fun_.extent(0); }
            inline SizeType n_quad_points() const override { return measure_.extent(1); }
            inline SizeType spatial_dimension() const override { return grad_.extent(3); }
            inline SizeType manifold_dimension() const override { return jacobian_.extent(1); }

            UTOPIA_INLINE_FUNCTION IntView &element_tags() { return element_tags_; }
            UTOPIA_INLINE_FUNCTION const IntView &element_tags() const { return element_tags_; }

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

            void print_measure() {
                printf("Measure:\n");
                auto n_quad_points = this->n_quad_points();
                auto n_cells = this->n_cells();

                // Avoids the capturing of 'this'
                auto measure = this->measure();

                Kokkos::parallel_for(
                    "FE::print_measure", n_cells, KOKKOS_LAMBDA(const int &cell) {
                        // if (cell != 0) return;

                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            printf("%g ", measure(cell, qp));
                        }

                        printf("\n");
                    });
            }

            void print_jacobian() {
                printf("Jacobian:\n");
                int n_cells = this->n_cells();
                int n_quad_points = this->n_quad_points();

                int spatial_dimension = this->spatial_dimension();
                int manifold_dimension = this->manifold_dimension();

                // Avoids the capturing of 'this'
                auto jacobian = this->jacobian();

                Kokkos::parallel_for(
                    "FE::print_jacobian", n_cells, KOKKOS_LAMBDA(const int &cell) {
                        // if (cell != 0) return;

                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            for (int r = 0; r < spatial_dimension; ++r) {
                                for (int c = 0; c < manifold_dimension; ++c) {
                                    printf("%g ", jacobian(cell, qp, r, c));
                                }
                                printf("\n");
                            }

                            printf("\n");
                        }

                        printf("\n");
                    });
            }

            void print_jacobian_inverse() {
                printf("Jacobian inverse:\n");
                int n_cells = this->n_cells();
                int n_quad_points = this->n_quad_points();

                int spatial_dimension = this->spatial_dimension();
                int manifold_dimension = this->manifold_dimension();

                // Avoids the capturing of 'this'
                auto jacobian_inv = this->jacobian_inverse();

                Kokkos::parallel_for(
                    "FE::print_jacobian_inverse", n_cells, KOKKOS_LAMBDA(const int &cell) {
                        // if (cell != 0) return;

                        printf("cell: %d\n", cell);

                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            for (int r = 0; r < manifold_dimension; ++r) {
                                for (int c = 0; c < spatial_dimension; ++c) {
                                    printf("%g ", jacobian_inv(cell, qp, r, c));
                                }
                                printf("\n");
                            }

                            printf("\n");
                        }

                        printf("\n");
                    });
            }

            void print_function() {
                printf("Function:\n");
                int n_quad_points = this->n_quad_points();
                int n_shape_functions = this->n_shape_functions();

                // Avoids the capturing of 'this'
                auto fun = this->fun();

                Kokkos::parallel_for(
                    "FE::print_function", 1, KOKKOS_LAMBDA(const int) {
                        for (int node = 0; node < n_shape_functions; ++node) {
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                printf("%g ", fun(node, qp));
                            }

                            printf("\n");
                        }
                    });
            }

            void print_gradient() {
                printf("Gradient:\n");

                int n_cells = this->n_cells();
                int n_quad_points = this->n_quad_points();

                int spatial_dimension = this->spatial_dimension();
                int n_shape_functions = this->n_shape_functions();

                // Avoids the capturing of 'this'
                auto gradient = this->grad();

                Kokkos::parallel_for(
                    "FE::print_gradient", n_cells, KOKKOS_LAMBDA(const int &cell) {
                        // if (cell != 0) return;

                        printf("cell: %d\n", cell);

                        for (int node = 0; node < n_shape_functions; ++node) {
                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                for (int d = 0; d < spatial_dimension; ++d) {
                                    printf("%g ", gradient(cell, node, qp, d));
                                }

                                printf("\n");
                            }

                            printf("\n");
                        }
                    });
            }

            inline bool has_element_tags() const { return element_tags_.size() > 0; }

        private:
            MeasureView measure_;
            FunctionView fun_;
            GradientView grad_;
            JacobianView jacobian_;
            JacobianInverseView jacobian_inverse_;

            IntView element_tags_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_FE_HPP
