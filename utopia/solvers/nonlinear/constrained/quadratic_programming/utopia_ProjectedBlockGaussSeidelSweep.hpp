#ifndef UTOPIA_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Trace.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Views.hpp"

#include "utopia_ProjectedGaussSeidelSweep.hpp"

#include <type_traits>
#include <vector>

namespace utopia {

    template <class Matrix, int BlockSize>
    class ProjectedBlockGaussSeidelSweep final : public ProjectedGaussSeidelSweep<Matrix> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        using ArrayView = utopia::ArrayView<Scalar>;
        using ConstArrayView = utopia::ArrayView<const Scalar>;
        using VectorView = utopia::VectorView<ArrayView>;
        using ConstVectorView = utopia::VectorView<ConstArrayView>;
        using Block = utopia::StaticMatrix<Scalar, BlockSize, BlockSize>;

        ProjectedBlockGaussSeidelSweep *clone() const override { return new ProjectedBlockGaussSeidelSweep(); }

        class BlockIdx {
        public:
            constexpr BlockIdx(const SizeType i_scalar, const SizeType j_scalar)
                : i(i_scalar / BlockSize),
                  sub_i(i_scalar - i * BlockSize),
                  j(j_scalar / BlockSize),
                  sub_j(j_scalar - j * BlockSize) {
                assert(sub_j < BlockSize);
                assert(sub_i < BlockSize);
            }

            SizeType i;
            SizeType sub_i;

            SizeType j;
            SizeType sub_j;

            inline constexpr bool is_diag() const { return i == j; }
        };

        void init_from_local_matrix(const Matrix &mat) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedBlockGaussSeidelSweep::init_from_local_matrix");
            const SizeType n_rows = mat.rows();
            const SizeType n_blocks = n_rows / BlockSize;

            if (n_blocks != SizeType(diag_.size())) {
                diag_.resize(n_blocks);
                inv_diag_.resize(n_blocks);
                row_ptr_.resize(n_rows + 1);
            }

            std::fill(row_ptr_.begin(), row_ptr_.end(), SizeType(0));

            SizeType n_off_diag_entries = 0;

            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                BlockIdx block(i, j);

                if (block.is_diag()) {
                    // store diagonal block (invert later)
                    diag_[block.i](block.sub_i, block.sub_j) = a_ij;

                    assert(a_ij != 0.0 || block.sub_i != block.sub_j);
                } else {
                    if (a_ij != 0.0) {
                        ++row_ptr_[i + 1];
                        ++n_off_diag_entries;
                    }
                }
            });

            for (SizeType b = 0; b < n_blocks; ++b) {
                assert(det(diag_[b]) > 0);
                inv_diag_[b] = inv(diag_[b]);
            }

            for (SizeType i = 0; i < n_rows; ++i) {
                row_ptr_[i + 1] += row_ptr_[i];
            }

            assert(n_off_diag_entries == row_ptr_.back());

            if (n_off_diag_entries != SizeType(values_.size())) {
                values_.resize(n_off_diag_entries);
                col_idx_.resize(n_off_diag_entries);
            }

            n_off_diag_entries = 0;
            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                BlockIdx block(i, j);

                if (!block.is_diag() && a_ij != 0.0) {
                    col_idx_[n_off_diag_entries] = j;
                    values_[n_off_diag_entries] = a_ij;
                    ++n_off_diag_entries;
                }
            });

            UTOPIA_TRACE_REGION_END("ProjectedBlockGaussSeidelSweep::init_from_local_matrix");
        }

        void update_from_local_matrix(const Matrix &local_diag_block) override {
            init_from_local_matrix(local_diag_block);
        }

        void apply(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedBlockGaussSeidelSweep::apply(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply();
            }

            UTOPIA_TRACE_REGION_END("ProjectedBlockGaussSeidelSweep::apply(...)");
        }

        void apply_unconstrained(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("ProjectedBlockGaussSeidelSweep::apply_unconstrained(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply_unconstrained();
            }

            UTOPIA_TRACE_REGION_END("ProjectedBlockGaussSeidelSweep::apply_unconstrained(...)");
        }

        ProjectedBlockGaussSeidelSweep() = default;

    private:
        void apply_unconstrained() {
            const SizeType n_blocks = diag_.size();

            StaticVector<Scalar, BlockSize> val, d_inv_val;

            auto f = [&](const SizeType &b) {
                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;

                    val[d] = this->r_[i];

                    const SizeType row_end = row_ptr_[i + 1];
                    for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val[d] -= values_[j] * this->c_[col_idx_[j]];
                    }
                }

                d_inv_val = inv_diag_[b] * val;

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    this->c_[i] = d_inv_val[d];
                }
            };

            for (SizeType b = 0; b < n_blocks; ++b) {
                f(b);
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType b = n_blocks - 1; b >= 0; --b) {
                    f(b);
                }
            }
        }

        void apply() {
            const SizeType n_blocks = diag_.size();

            StaticVector<Scalar, BlockSize> val, d_inv_val;
            StaticMatrix<Scalar, BlockSize, BlockSize> d_inv, d_temp;

            auto f = [&](const SizeType &b) {
                d_temp.copy(diag_[b]);

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;

                    val[d] = this->r_[i];

                    const SizeType row_end = row_ptr_[i + 1];
                    for (SizeType j = row_ptr_[i]; j < row_end; ++j) {
                        val[d] -= values_[j] * this->c_[col_idx_[j]];
                    }
                }

                d_inv_val = inv_diag_[b] * val;

                bool constrained = false;
                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    if (d_inv_val[d] > this->lb_[i] && d_inv_val[d] < this->ub_[i]) continue;

                    constrained = true;

                    val[d] = device::max(this->lb_[i], device::min(d_inv_val[d], this->ub_[i]));

                    for (SizeType d_j = 0; d_j < BlockSize; ++d_j) {
                        d_temp(d, d_j) = d_j == d;
                    }
                }

                if (constrained) {
                    d_inv = inv(d_temp);
                    d_inv_val = d_inv * val;
                }

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    this->c_[i] = d_inv_val[d];
                }
            };

            for (SizeType b = 0; b < n_blocks; ++b) {
                f(b);
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType b = n_blocks - 1; b >= 0; --b) {
                    f(b);
                }
            }
        }

        std::vector<Block> diag_;
        std::vector<Block> inv_diag_;
        std::vector<Scalar> values_;
        std::vector<SizeType> row_ptr_;
        std::vector<SizeType> col_idx_;
    };
}  // namespace utopia

#endif  // UTOPIA_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
