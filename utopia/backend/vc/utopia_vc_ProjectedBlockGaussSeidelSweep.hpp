#ifndef UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_VectorView.hpp"
#include "utopia_Views.hpp"

#include "utopia_CRSMatrix.hpp"
#include "utopia_CRSToBlockCRS.hpp"
#include "utopia_ProjectedBlockGaussSeidelSweep.hpp"

#include "utopia_petsc_ILUDecompose.hpp"

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
        using ArrayViewT = utopia::ArrayView<Scalar, DYNAMIC_SIZE, DYNAMIC_SIZE>;
        using BlockView = utopia::TensorView<ArrayViewT, 2>;

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
            crs_block_matrix<BlockSize>(mat, block_crs);

            auto n_blocks = block_crs.rows();

            inv_diag_.resize(n_blocks);
            diag_.resize(n_blocks);

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                auto row = block_crs.row(block_i);

                for (SizeType k = 0; k < row.n_blocks(); ++k) {
                    if (row.colidx(k) == block_i) {
                        d.raw_type().set_data(row.block(k));
                        diag_[block_i].copy(d);
                        inv_diag_[block_i] = inv(diag_[block_i]);
                        break;
                    }
                }
            }

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
            auto &row_ptr = block_crs.row_ptr();
            auto &colidx = block_crs.colidx();

            const SizeType n_blocks = row_ptr.size() - 1;

            SIMDType r_simd, c_simd, LRc_simd, mat_simd;

            BlockView aij;
            aij.raw_type().set_size(BlockSize, BlockSize);

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b_offset + d;
                    r_simd[d] = this->r_[i];
                }

                for (SizeType j = row_ptr[block_i]; j < row_end; ++j) {
                    if (colidx[j] == block_i) continue;

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        c_simd[d] = this->c_[colidx[j] * BlockSize + d];
                    }

                    aij.raw_type().set_data(block_crs.block(j));

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&(aij(d, 0)), Vc::Unaligned);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[block_i];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), Vc::Unaligned);
                    SizeType i = b_offset + d;
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
            auto &row_ptr = block_crs.row_ptr();
            auto &colidx = block_crs.colidx();

            const SizeType n_blocks = row_ptr.size() - 1;

            SmallMatrix d_temp, d_inv_temp;
            SIMDType r_simd, c_simd, LRc_simd, mat_simd;

            BlockView aij;
            aij.raw_type().set_size(BlockSize, BlockSize);

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b_offset + d;
                    r_simd[d] = this->r_[i];
                }

                for (SizeType j = row_ptr[block_i]; j < row_end; ++j) {
                    if (colidx[j] == block_i) continue;

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        c_simd[d] = this->c_[colidx[j] * BlockSize + d];
                    }

                    aij.raw_type().set_data(block_crs.block(j));

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&(aij(d, 0)), Vc::Unaligned);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[block_i];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), Vc::Unaligned);
                    SizeType i = b_offset + d;
                    this->c_[i] = (mat_simd * r_simd).sum();
                }

                d_temp.copy(diag_[block_i]);

                bool constrained = false;
                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b_offset + d;
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
                        SizeType i = b_offset + d;
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
        CRSMatrix<std::vector<Scalar>, std::vector<SizeType>, BlockSize> block_crs;
    };
}  // namespace utopia

#endif  // UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
