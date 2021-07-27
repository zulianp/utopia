#ifndef UTOPIA_INTREPID2_RESIDUAL_HPP
#define UTOPIA_INTREPID2_RESIDUAL_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_intrepid2_Base.hpp"

namespace utopia {

    namespace intrepid2 {

        template <typename T, typename... Args>
        void residual(const View<T, Args...> &cell_mat,
                      const View<T, Args...> &cell_vec,
                      View<T, Args...> &cell_result,
                      AssemblyMode mode) {
            using ExecutionSpace = typename View<T, Args...>::execution_space;
            using CellTestRange = Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>;
            using SizeType = ::Intrepid2::ordinal_type;

            const SizeType n_cells = cell_mat.extent(0);
            const SizeType n_nodes_x_cell = cell_mat.extent(1);
            const SizeType n = cell_mat.extent(2);

            switch (mode) {
                case ADD_MODE: {
                    ::Kokkos::parallel_for(
                        "residual",
                        CellTestRange({0, 0}, {n_cells, n_nodes_x_cell}),
                        UTOPIA_LAMBDA(const SizeType cell, const SizeType i) {
                            for (SizeType j = 0; j < n; ++j) {
                                cell_result(cell, i) += cell_mat(cell, i, j) * cell_vec(cell, j);
                            }
                        });
                    break;
                }
                case SUBTRACT_MODE: {
                    ::Kokkos::parallel_for(
                        "residual",
                        CellTestRange({0, 0}, {n_cells, n_nodes_x_cell}),
                        UTOPIA_LAMBDA(const SizeType cell, const SizeType i) {
                            for (SizeType j = 0; j < n; ++j) {
                                cell_result(cell, i) -= cell_mat(cell, i, j) * cell_vec(cell, j);
                            }
                        });

                    break;
                }
                case OVERWRITE_MODE: {
                    ::Kokkos::parallel_for(
                        "residual",
                        CellTestRange({0, 0}, {n_cells, n_nodes_x_cell}),
                        UTOPIA_LAMBDA(const SizeType cell, const SizeType i) {
                            cell_result(cell, i) = 0;
                            for (SizeType j = 0; j < n; ++j) {
                                cell_result(cell, i) -= cell_mat(cell, i, j) * cell_vec(cell, j);
                            }
                        });
                    break;
                }
                default: {
                    assert(false);
                    Utopia::Abort("Mode does not exist!");
                }
            }
        }
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_RESIDUAL_HPP