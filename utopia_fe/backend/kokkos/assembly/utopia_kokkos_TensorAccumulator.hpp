#ifndef UTOPIA_KOKKOS_TENSOR_ACCUMULATOR_HPP
#define UTOPIA_KOKKOS_TENSOR_ACCUMULATOR_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

#include "utopia_kokkos_Commons.hpp"
#include "utopia_kokkos_Traits.hpp"

namespace utopia {
    namespace kokkos {

        template <class View_>
        class TensorAccumulator : public Describable {
        public:
            using View = View_;
            using Scalar = typename utopia::Traits<View>::Scalar;
            // using ExecutionSpace = typename Traits<View>::ExecutionSpace;
            using ExecutionSpace = typename View::execution_space;

            View &data() { return data_; }

            template <class FE>
            bool is_compatible(const FE &fe) const {
                return fe.n_cells() <= data_.extent(0);
            }

            template <class FE>
            void init_scalar(const FE &fe, int n_vars) {
                data_ = View("scalars", fe.n_cells(), n_vars);
            }

            template <class FE>
            void init_matrix(const FE &fe, int n_vars) {
                const int n_shape_functions = fe.n_shape_functions();
                const int n_dofs = n_shape_functions * n_vars;
                data_ = View("matrices", fe.n_cells(), n_dofs, n_dofs);
            }

            template <class FE>
            void init_vector(FE &fe, int n_vars) {
                const int n_shape_functions = fe.n_shape_functions();
                const int n_dofs = n_shape_functions * n_vars;
                data_ = View("vectors", fe.n_cells(), n_dofs);
            }

            inline AssemblyMode mode() const { return mode_; }
            inline void set_mode(AssemblyMode mode) const { mode_ = mode; }

            void prepare() {
                if (mode_ == OVERWRITE_MODE) {
                    zero();
                }
            }

            Scalar sum() const {
                auto data = data_.data();
                auto n_elements = data_.size();

                Scalar ret = 0.0;
                Kokkos::parallel_reduce(
                    Kokkos::RangePolicy<ExecutionSpace>(0, n_elements),
                    UTOPIA_LAMBDA(int i, Scalar &val) { val += data[i]; },
                    ret);

                return ret;
            }

            void zero() { utopia::kokkos::fill(data_, 0.0); }

            void describe(std::ostream &os) const override {
                const SizeType n_cells = data_.extent(0);

                if (data_.rank() == 2) {
                    const int n_dofs_i = data_.extent(1);

                    for (SizeType c = 0; c < n_cells; ++c) {
                        os << c << ")\n";
                        for (SizeType i = 0; i < n_dofs_i; ++i) {
                            os << data_(c, i) << " ";

                            os << '\n';
                        }

                        os << '\n';
                    }

                } else if (data_.rank() == 3) {
                    const int n_dofs_i = data_.extent(1);
                    const int n_dofs_j = data_.extent(2);

                    for (SizeType c = 0; c < n_cells; ++c) {
                        os << c << ")\n";
                        for (SizeType i = 0; i < n_dofs_i; ++i) {
                            for (SizeType j = 0; j < n_dofs_j; ++j) {
                                os << data_(c, i, j) << " ";
                            }

                            os << '\n';
                        }

                        os << '\n';
                    }
                }
            }

        private:
            View data_;
            AssemblyMode mode_{ADD_MODE};
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TENSOR_ACCUMULATOR_HPP