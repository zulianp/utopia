#if 1
#include "utopia_vc_ProjectedBlockGaussSeidelSweep_crahes_on_daint.hpp"
#else

#ifndef UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
#define UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_PETSC

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
        static const int BlockSize_2 = BlockSize * BlockSize;

        using SmallMatrix = utopia::StaticMatrix<Scalar, BlockSize, BlockSize>;
        using Block = SmallMatrix;  // utopia::StaticVector<SIMDType, BlockSize>;
        using ArrayViewT = utopia::ArrayView<Scalar, DYNAMIC_SIZE, DYNAMIC_SIZE>;
        using BlockView = utopia::TensorView<ArrayViewT, 2>;

        static_assert(BlockSize == 4, "Must fix this class to work with other types than avx2 with doubles");

        VcProjectedBlockGaussSeidelSweep *clone() const override {
            auto ptr = utopia::make_unique<VcProjectedBlockGaussSeidelSweep>();
            ptr->on_update_keep_nnz_pattern_ = on_update_keep_nnz_pattern_;
            return ptr.release();
        }

        void read(Input &in) override { in.get("on_update_keep_nnz_pattern", on_update_keep_nnz_pattern_); }

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
            if (diag_.empty() || !on_update_keep_nnz_pattern_) {
                UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweep::init_crs");
                crs_block_matrix_split_diag<BlockSize>(mat, block_crs, diag_);
                inv_diag_.resize(block_crs.rows());
                UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweep::init_crs");
            } else {
                UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweep::update_crs");
                crs_block_matrix_update<BlockSize>(mat, block_crs, diag_);
                UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweep::update_crs");
            }

            auto n_blocks = block_crs.rows();

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                d.raw_type().set_data(&diag_[block_i * BlockSize_2]);
                assert(std::abs(det(d)) > 0);
                inv_diag_[block_i] = inv(d);
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

            static const auto alignment = Vc::Unaligned;
            // static const auto alignment = Vc::Aligned;

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                r_simd.load(&this->r_[b_offset], alignment);

                for (SizeType j = row_ptr[block_i]; j < row_end; ++j) {
                    c_simd.load(&this->c_[colidx[j] * BlockSize], alignment);

                    auto *aij = block_crs.block(j);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&aij[d * BlockSize], alignment);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[block_i];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), alignment);
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

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            static const auto alignment = Vc::Unaligned;
            // static const auto alignment = Vc::Aligned;
            bool coord_is_constrained[BlockSize];

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                r_simd.load(&this->r_[b_offset], alignment);

                for (SizeType j = row_ptr[block_i]; j < row_end; ++j) {
                    c_simd.load(&this->c_[colidx[j] * BlockSize], alignment);

                    auto *aij = block_crs.block(j);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&aij[d * BlockSize], alignment);
                        LRc_simd[d] = (mat_simd * c_simd).sum();
                    }

                    r_simd -= LRc_simd;
                }

                const auto &D_inv = inv_diag_[block_i];
                for (SizeType d = 0; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), alignment);
                    SizeType i = b_offset + d;
                    this->c_[i] = (mat_simd * r_simd).sum();
                }

                bool constrained = false;
                for (SizeType d = 0; d < BlockSize; ++d) {
                    SizeType i = b_offset + d;
                    auto c_i = this->c_[i];

                    // Yes it is an assignment...
                    if ((coord_is_constrained[d] = !(c_i > this->lb_[i] && c_i < this->ub_[i]))) {
                        constrained = true;
                        r_simd[d] = device::max(this->lb_[i], device::min(c_i, this->ub_[i]));
                    }
                }

                if (constrained) {
                    d.raw_type().set_data(&diag_[block_i * BlockSize_2]);
                    d_temp.copy(d);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        if (coord_is_constrained[d]) {
                            for (SizeType d_j = 0; d_j < BlockSize; ++d_j) {
                                d_temp(d, d_j) = d_j == d;
                            }
                        }
                    }

                    d_inv_temp = inv(d_temp);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        SizeType i = b_offset + d;
                        mat_simd.load(&(d_inv_temp(d, 0)), alignment);
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

        void on_update_keep_nnz_pattern(const bool val) { on_update_keep_nnz_pattern_ = val; }

    private:
        std::vector<Scalar> diag_;
        std::vector<Block> inv_diag_;
        CRSMatrix<std::vector<Scalar>, std::vector<SizeType>, BlockSize> block_crs;
        bool on_update_keep_nnz_pattern_{false};
    };
}  // namespace utopia

#endif
#endif  // UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_HPP
#endif
