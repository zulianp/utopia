#ifndef UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Views.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Algorithms.hpp"

#include <vector>
#include <type_traits>

namespace utopia {

    template<typename Scalar, typename SizeType>
    class ProjectedGaussSeidelSweep {
    public:
        using ArrayView  = utopia::ArrayView<Scalar>;
        using ConstArrayView  = utopia::ArrayView<const Scalar>;
        using VectorView = utopia::VectorView<ArrayView>;
        using ConstVectorView = utopia::VectorView<ConstArrayView>;

        template<class Derived>
        void init_from_local_matrix(const Tensor<Derived, 2> &local_diag_block)
        {
            const auto &mat = local_diag_block.derived();

            const SizeType n_rows = mat.rows();

            if(n_rows != SizeType(d_inv_.size())) {
                d_inv_.resize(n_rows);
                row_ptr_.resize(n_rows + 1);
            }

            std::fill(row_ptr_.begin(), row_ptr_.end(), SizeType(0));

            SizeType n_off_diag_entries = 0;

            mat.read(
                [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                    if(i == j) {
                        d_inv_[i] = (device::abs(a_ij) > 0.0)? (1/a_ij) : Scalar(0.0);
                    } else {
                        ++row_ptr_[i+1];
                        ++n_off_diag_entries;
                    }
            });


            for(SizeType i = 0; i < n_rows; ++i) {
                row_ptr_[i+1] += row_ptr_[i];
            }

            assert(n_off_diag_entries == row_ptr_.back());


            if(n_off_diag_entries != SizeType(values_.size())) {
                values_.resize(n_off_diag_entries);
                col_idx_.resize(n_off_diag_entries);
            }


            n_off_diag_entries = 0;
            mat.read(
            [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                    if(i != j) {
                        col_idx_[n_off_diag_entries] = j;
                        values_[n_off_diag_entries] = a_ij;
                        ++n_off_diag_entries;
                    }
            });
        }


        template<class Derived>
        void update_from_local_matrix(const Tensor<Derived, 2> &local_diag_block)
        {
            init_from_local_matrix(local_diag_block);

            // const auto &mat = local_diag_block.derived();
            // const SizeType n_rows = local_diag_block.rows();

            // d_inv_.resize(n_rows);
            // row_ptr_.resize(n_rows + 1);

            // row_ptr_[0] = 0;
            // SizeType n_off_diag_entries = 0;

            // local_diag_block.read(
            //     [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
            //         if(i == j) {
            //             d_inv_[i] = (device::abs(a_ij) > 0.0)? (1/a_ij) : Scalar(0.0);
            //         } else {
            //             ++row_ptr_[i+1];
            //             ++n_off_diag_entries;
            //         }
            // });

            // values_.resize(n_off_diag_entries);
            // col_idx_.resize(n_off_diag_entries);

            // SizeType n_off_diag_entries = 0;
            // mat.read(
            // [&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
            //         if(i == j) {
            //             d_inv_[i] = (device::abs(a_ij) > 0.0)? (1/a_ij) : Scalar(0.0);
            //         } else {
                        // col_idx_[n_off_diag_entries] = j;
            //             assert(n_off_diag_entries < SizeType(values_.size()));
            //             values_[n_off_diag_entries] = a_ij;
            //             ++n_off_diag_entries;
            //         }
            // });
        }

        void apply(const SizeType &times)
        {
            for(SizeType t = 0; t < times; ++t) {
                apply();
            }
        }

        void apply_unconstrained(const SizeType &times)
        {
            for(SizeType t = 0; t < times; ++t) {
                apply_unconstrained();
            }
        }

        void apply_unconstrained()
        {
            const SizeType n_rows = d_inv_.size();

            for(SizeType i = 0; i < n_rows; ++i) {
                Scalar val = r_[i];

                const SizeType row_end = row_ptr_[i+1];
                for(SizeType j = row_ptr_[i]; j < row_end; ++j) {
                    val -= values_[j] * c_[col_idx_[j]];
                }

                c_[i] = d_inv_[i] * val;
            }

            if(symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for(SizeType i = n_rows-1; i >= 0; --i) {
                    Scalar val = r_[i];

                    const SizeType row_end = row_ptr_[i+1];
                    for(SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val -= values_[j] * c_[col_idx_[j]];
                    }

                    c_[i] = d_inv_[i] * val;
                }
            }
        }

        void apply()
        {
            const SizeType n_rows = d_inv_.size();

            for(SizeType i = 0; i < n_rows; ++i) {
                Scalar val = r_[i];

                const SizeType row_end = row_ptr_[i+1];
                for(SizeType j = row_ptr_[i]; j < row_end; ++j) {
                    val -= values_[j] * c_[col_idx_[j]];
                }

                c_[i] = device::max(
                    lb_[i],
                    device::min(
                        d_inv_[i] * val,
                        ub_[i]
                    )
                );
            }

            if(symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for(SizeType i = n_rows-1; i >= 0; --i) {
                    Scalar val = r_[i];

                    const SizeType row_end = row_ptr_[i+1];
                    for(SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val -= values_[j] * c_[col_idx_[j]];
                    }

                    c_[i] = device::max(
                        lb_[i],
                        device::min(
                            d_inv_[i] * val,
                            ub_[i]
                        )
                    );
                }
            }
        }

        void set_residual_view(const ConstVectorView &r)
        {
            //FIXME
            r_.raw_type().set(&r[0], r.size());
        }

        void set_correction_view(VectorView &c)
        {
            //FIXME
            c_.raw_type().set(&c[0], c.size());
        }

        void set_correction_view(VectorView &&c)
        {
            //FIXME
            c_.raw_type().set(&c[0], c.size());
        }

        void set_bounds(const ConstVectorView &lb, const ConstVectorView &ub)
        {
            //FIXME
            lb_.raw_type().set(&lb[0], lb.size());
            ub_.raw_type().set(&ub[0], ub.size());
        }

        void symmetric(const bool symmetric)
        {
            symmetric_ = symmetric_;
        }

        void l1(const bool l1)
        {
            l1_ = l1;
        }

        ProjectedGaussSeidelSweep()
        : symmetric_(true)
        {}

    private:
        VectorView c_;
        ConstVectorView r_, lb_, ub_;

        std::vector<Scalar> values_;
        std::vector<Scalar> d_inv_;
        std::vector<SizeType> row_ptr_;
        std::vector<SizeType> col_idx_;

        bool symmetric_;
        bool l1_;
    };
}


#endif //UTOPIA_PROJECTED_GAUSS_SEIDEL_SWEEP_HPP
