#ifndef UTOPIA_KOKKOS_FE_BASE_HPP
#define UTOPIA_KOKKOS_FE_BASE_HPP

#include <Kokkos_DynRankView.hpp>

namespace utopia {
    namespace kokkos {

        template <typename T>
        using DefaultView = ::Kokkos::DynRankView<T, ::Kokkos::DefaultExecutionSpace>;

        template <typename ExecutionSpace_, typename SizeType>
        class FEBase {
        public:
            using ExecutionSpace = ExecutionSpace_;
            using CellTestTrialRange = Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>;
            using CellTestRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;
            using CellRange = Kokkos::RangePolicy<ExecutionSpace>;
            using CellQPRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;

            UTOPIA_FUNCTION ~FEBase() = default;

            virtual SizeType n_cells() const = 0;
            virtual SizeType n_shape_functions() const = 0;
            virtual SizeType n_quad_points() const = 0;
            virtual SizeType spatial_dimension() const = 0;
            virtual SizeType manifold_dimension() const = 0;

            virtual bool is_shell() const { return manifold_dimension() < spatial_dimension(); }

            inline CellQPRange cell_qp_range() const {
                SizeType n_cells = this->n_cells();
                SizeType n_quad_points = this->n_quad_points();
                return CellQPRange({0, 0}, {n_cells, n_quad_points});
            }

            static inline Kokkos::RangePolicy<ExecutionSpace> range(const SizeType begin, const SizeType end) {
                return Kokkos::RangePolicy<ExecutionSpace>(begin, end);
            }

            inline CellRange cell_range() const {
                SizeType n_cells = this->n_cells();
                return CellRange(0, n_cells);
            }

            inline CellTestTrialRange cell_test_trial_range() const {
                SizeType n_cells = this->n_cells();
                SizeType n_shape_functions = this->n_shape_functions();

                return CellTestTrialRange({0, 0, 0}, {n_cells, n_shape_functions, n_shape_functions});
            }

            inline CellTestRange cell_test_range() const {
                SizeType n_cells = this->n_cells();
                SizeType n_shape_functions = this->n_shape_functions();
                return CellTestRange({0, 0}, {n_cells, n_shape_functions});
            }
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FE_BASE_HPP
