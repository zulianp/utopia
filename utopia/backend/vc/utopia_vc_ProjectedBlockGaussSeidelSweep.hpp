#ifndef UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Views.hpp"

#include "utopia_ProjectedBlockGaussSeidelSweep.hpp"

#include <map>
#include <type_traits>
#include <vector>

namespace utopia {

    template <class Matrix>
    using MatrixSIMDType = Vc::Vector<typename Traits<Matrix>::Scalar>;

    // Slower than original non-vectorized version
    template <class Matrix>
    class VcProjectedBlockGaussSeidelSweep final : public ProjectedGaussSeidelSweep<Matrix> {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        using ArrayView = utopia::ArrayView<Scalar>;
        using ConstArrayView = utopia::ArrayView<const Scalar>;
        using VectorView = utopia::VectorView<ArrayView>;
        using ConstVectorView = utopia::VectorView<ConstArrayView>;

        using SIMDType = MatrixSIMDType<Matrix>;
        static const int BlockSize = static_cast<int>(SIMDType::Size);

        using SmallMatrix = utopia::StaticMatrix<Scalar, BlockSize, BlockSize>;
        using Block = SmallMatrix;  // utopia::StaticVector<SIMDType, BlockSize>;

        static_assert(BlockSize == 4, "Must fix this class to work with other types than avx2 with doubles");

        VcProjectedBlockGaussSeidelSweep *clone() const override { return new VcProjectedBlockGaussSeidelSweep(); }

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
            inline constexpr bool is_first() const { return sub_i == 0 && sub_j == 0; }
        };

        void init_from_local_matrix(const Matrix &mat) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweep::init_from_local_matrix");

            const SizeType n_rows = mat.rows();
            const SizeType n_blocks = n_rows / BlockSize;

            assert(n_blocks * BlockSize == n_rows && "Rows must be a multiple of BlockSize");

            if (n_blocks != SizeType(diag_.size())) {
                diag_.resize(n_blocks);
                inv_diag_.resize(n_blocks);
                row_ptr_.resize(n_blocks + 1);
            }

            std::fill(row_ptr_.begin(), row_ptr_.end(), SizeType(0));

            SizeType n_off_diag_entries = 0;

            // FIXME use unorderd map which should be faster
            // std::unordered_map<std::pair<SizeType, SizeType>, SizeType> block_counted;
            std::map<std::pair<SizeType, SizeType>, SizeType> block_counted;

            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                BlockIdx block(i, j);

                assert(block.i < SizeType(diag_.size()));

                if (block.is_diag()) {
                    diag_[block.i](block.sub_i, block.sub_j) = a_ij;

                    assert(a_ij != 0.0 || block.sub_i != block.sub_j);
                } else {
                    auto ij = std::make_pair(std::pair<SizeType, SizeType>(block.i, block.j), n_off_diag_entries);
                    auto ret = block_counted.insert(ij);

                    if (ret.second) {
                        ++row_ptr_[block.i + 1];
                        ++n_off_diag_entries;
                    }
                }
            });

            for (SizeType b = 0; b < n_blocks; ++b) {
                assert(det(diag_[b]) > 0);
                inv_diag_[b] = inv(diag_[b]);
            }

            for (SizeType i = 0; i < n_blocks; ++i) {
                row_ptr_[i + 1] += row_ptr_[i];
            }

            assert(n_off_diag_entries == row_ptr_.back());

            if (n_off_diag_entries != SizeType(values_.size())) {
                values_.resize(n_off_diag_entries);
                col_idx_.resize(n_off_diag_entries);
            }

            for (auto &A : values_) {
                A.set(Scalar(0.0));
            }

            mat.read([&](const SizeType &i, const SizeType &j, const Scalar &a_ij) {
                BlockIdx block(i, j);

                if (!block.is_diag()) {
                    auto it = block_counted.find(std::pair<SizeType, SizeType>(block.i, block.j));

                    assert(it->second < SizeType(col_idx_.size()));
                    col_idx_[it->second] = block.j;
                    values_[it->second](block.sub_i, block.sub_j) = a_ij;
                }
            });

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweep::init_from_local_matrix");
        }

        void update_from_local_matrix(const Matrix &local_diag_block) override {
            init_from_local_matrix(local_diag_block);
        }

        void apply(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweep::apply(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply();
            }

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweep::apply(...)");
        }

        void apply_unconstrained(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweep::apply_unconstrained(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply_unconstrained();
            }

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweep::apply_unconstrained(...)");
        }

        VcProjectedBlockGaussSeidelSweep() = default;

    private:
        void apply_unconstrained() {
            const SizeType n_blocks = row_ptr_.size() - 1;

            SIMDType r_simd, c_simd, LRc_simd, mat_simd;

            auto f = [&](const SizeType &b) {
                const SizeType row_end = row_ptr_[b + 1];

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    r_simd[d] = this->r_[i];
                }

                for (SizeType j = row_ptr_[b]; j < row_end; ++j) {
                    for (SizeType d = 0; d < BlockSize; ++d) {
                        c_simd[d] = this->c_[col_idx_[j] * BlockSize + d];
                    }

                    const auto &LR = values_[j];
                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&(LR(d, 0)), Vc::Unaligned);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[b];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), Vc::Unaligned);
                    SizeType i = b * BlockSize + d;
                    this->c_[i] = (mat_simd * r_simd).sum();
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
            const SizeType n_blocks = row_ptr_.size() - 1;

            SmallMatrix d_temp, d_inv_temp;
            SIMDType r_simd, c_simd, LRc_simd, mat_simd;

            auto f = [&](const SizeType &b) {
                const SizeType row_end = row_ptr_[b + 1];

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    r_simd[d] = this->r_[i];
                }

                for (SizeType j = row_ptr_[b]; j < row_end; ++j) {
                    for (SizeType d = 0; d < BlockSize; ++d) {
                        c_simd[d] = this->c_[col_idx_[j] * BlockSize + d];
                    }

                    const auto &LR = values_[j];
                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&(LR(d, 0)), Vc::Unaligned);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[b];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), Vc::Unaligned);
                    SizeType i = b * BlockSize + d;
                    this->c_[i] = (mat_simd * r_simd).sum();
                }

                d_temp.copy(diag_[b]);

                bool constrained = false;
                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b * BlockSize + d;
                    auto c_i = this->c_[i];
                    if (c_i > this->lb_[i] && c_i < this->ub_[i]) continue;

                    constrained = true;

                    r_simd[d] = device::max(this->lb_[i], device::min(c_i, this->ub_[i]));

                    for (SizeType d_j = 0; d_j < BlockSize; ++d_j) {
                        d_temp(d, d_j) = d_j == d;
                    }
                }

                if (constrained) {
                    d_inv_temp = inv(d_temp);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        SizeType i = b * BlockSize + d;
                        mat_simd.load(&(d_inv_temp(d, 0)), Vc::Unaligned);
                        this->c_[i] = (mat_simd * r_simd).sum();
                    }
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
        std::vector<Block> values_;
        std::vector<SizeType> row_ptr_;
        std::vector<SizeType> col_idx_;
    };
}  // namespace utopia

#endif  // UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
