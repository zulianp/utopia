#ifndef UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_TRANSPOSED_HPP
#define UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_TRANSPOSED_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_ENABLE_PETSC

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
    class VcProjectedBlockGaussSeidelSweepTransposed final : public ProjectedGaussSeidelSweep<Matrix> {
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

        VcProjectedBlockGaussSeidelSweepTransposed *clone() const override {
            auto ptr = utopia::make_unique<VcProjectedBlockGaussSeidelSweepTransposed>();
            ptr->on_update_keep_nnz_pattern_ = on_update_keep_nnz_pattern_;
            return ptr.release();
        }

        void read(Input &in) override { in.get("on_update_keep_nnz_pattern", on_update_keep_nnz_pattern_); }

        void init_from_local_matrix(const Matrix &mat) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweepTransposed::init_from_local_matrix");
            if (diag_.empty() || !on_update_keep_nnz_pattern_) {
                UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweepTransposed::init_crs");
                crs_block_matrix_split_diag<BlockSize>(mat, block_crs, diag_);
                inv_diag_.resize(block_crs.rows());
                UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweepTransposed::init_crs");
            } else {
                UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweepTransposed::update_crs");
                crs_block_matrix_update<BlockSize>(mat, block_crs, diag_);
                UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweepTransposed::update_crs");
            }

            auto n_block_rows = block_crs.rows();
            auto n_blocks = block_crs.n_blocks();

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            for (SizeType k = 0; k < n_blocks; ++k) {
                d.raw_type().set_data(block_crs.block(k));

                for (int d1 = 0; d1 < BlockSize; ++d1) {
                    for (int d2 = d1 + 1; d2 < BlockSize; ++d2) {
                        std::swap(d(d1, d2), d(d2, d1));
                    }
                }
            }

            for (SizeType block_i = 0; block_i < n_block_rows; ++block_i) {
                d.raw_type().set_data(&diag_[block_i * BlockSize_2]);

                for (int d1 = 0; d1 < BlockSize; ++d1) {
                    for (int d2 = d1 + 1; d2 < BlockSize; ++d2) {
                        std::swap(d(d1, d2), d(d2, d1));
                    }
                }

                assert(std::abs(det(d)) > 0);
                inv_diag_[block_i] = inv(d);
            }

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweepTransposed::init_from_local_matrix");
        }

        void update_from_local_matrix(const Matrix &local_diag_block) override {
            init_from_local_matrix(local_diag_block);
        }

        void apply(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweepTransposed::apply(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply();
            }

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweepTransposed::apply(...)");
        }

        void apply_unconstrained(const SizeType &times) override {
            UTOPIA_TRACE_REGION_BEGIN("VcProjectedBlockGaussSeidelSweepTransposed::apply_unconstrained(...)");

            for (SizeType t = 0; t < times; ++t) {
                apply_unconstrained();
            }

            UTOPIA_TRACE_REGION_END("VcProjectedBlockGaussSeidelSweepTransposed::apply_unconstrained(...)");
        }

        VcProjectedBlockGaussSeidelSweepTransposed() = default;

    private:
        void apply_unconstrained() {
            auto &row_ptr = block_crs.row_ptr();
            auto &colidx = block_crs.colidx();

            const SizeType n_block_rows = row_ptr.size() - 1;

            SIMDType r_simd, c_simd, mat_simd;
            Scalar r_v[BlockSize] = {0, 0, 0, 0};

            SIMDType temp_c[BlockSize];

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            static const auto alignment = Vc::Unaligned;
            // static const auto alignment = Vc::Aligned;

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                r_simd.load(&this->r_[b_offset], alignment);

                for (SizeType k = row_ptr[block_i]; k < row_end; ++k) {
                    auto *aij = block_crs.block(k);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        temp_c[d] = this->c_[colidx[k] * BlockSize + d];
                    }

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&aij[d * BlockSize], alignment);
                        r_simd -= (temp_c[d] * mat_simd);
                    }
                }

                const auto &D_inv = inv_diag_[block_i];

                ////////////////////////////////////////////////

                mat_simd.load(&(D_inv(0, 0)), alignment);
                r_simd.store(r_v, alignment);

                c_simd = mat_simd * r_v[0];

                for (SizeType d = 1; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), alignment);
                    c_simd += mat_simd * r_v[d];
                }

                c_simd.store(&this->c_[b_offset], alignment);
            };

            for (SizeType b = 0; b < n_block_rows; ++b) {
                f(b);
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType b = n_block_rows - 1; b >= 0; --b) {
                    f(b);
                }
            }
        }

        void apply() {
            auto &row_ptr = block_crs.row_ptr();
            auto &colidx = block_crs.colidx();

            const SizeType n_block_rows = row_ptr.size() - 1;

            SmallMatrix d_temp, d_inv_temp;
            SIMDType r_simd, c_simd, mat_simd, l_simd, u_simd;
            Scalar r_v[BlockSize] = {0, 0, 0, 0};

            BlockView d;
            d.raw_type().set_size(BlockSize, BlockSize);

            static const auto alignment = Vc::Unaligned;
            // static const auto alignment = Vc::Aligned;

            auto f = [&](const SizeType &block_i) {
                const SizeType row_end = row_ptr[block_i + 1];
                const SizeType b_offset = block_i * BlockSize;

                r_simd.load(&this->r_[b_offset], alignment);

                for (SizeType k = row_ptr[block_i]; k < row_end; ++k) {
                    auto *aij = block_crs.block(k);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        mat_simd.load(&aij[d * BlockSize], alignment);
                        r_simd -= (this->c_[colidx[k] * BlockSize + d] * mat_simd);
                    }
                }

                const auto &D_inv = inv_diag_[block_i];

                ////////////////////////////////////////////////

                mat_simd.load(&(D_inv(0, 0)), alignment);
                r_simd.store(r_v, alignment);

                c_simd = mat_simd * r_v[0];

                for (SizeType d = 1; d < BlockSize; ++d) {
                    mat_simd.load(&(D_inv(d, 0)), alignment);
                    c_simd += mat_simd * r_v[d];
                }

                ////////////////////////////////////////////////

                l_simd.load(&this->lb_[b_offset], alignment);
                u_simd.load(&this->ub_[b_offset], alignment);

                auto mask_l = c_simd <= l_simd;
                auto mask_u = c_simd >= u_simd;
                auto mask_r = mask_l | mask_u;

                if (mask_r.isNotEmpty()) {
                    l_simd.setZeroInverted(mask_l);
                    u_simd.setZeroInverted(mask_u);
                    r_simd.setZero(mask_r);
                    r_simd += l_simd + u_simd;

                    d.raw_type().set_data(&diag_[block_i * BlockSize_2]);
                    d_temp.copy(d);

                    for (SizeType d = 0; d < BlockSize; ++d) {
                        if (mask_r[d]) {
                            for (SizeType d_j = 0; d_j < BlockSize; ++d_j) {
                                d_temp(d_j, d) = d_j == d;
                            }
                        }
                    }

                    d_inv_temp = inv(d_temp);

                    ////////////////////////////////////////////////
                    mat_simd.load(&(d_inv_temp(0, 0)), alignment);
                    r_simd.store(r_v, alignment);

                    c_simd = mat_simd * r_v[0];

                    for (SizeType d = 1; d < BlockSize; ++d) {
                        mat_simd.load(&(d_inv_temp(d, 0)), alignment);
                        c_simd += mat_simd * r_v[d];
                    }

                    ////////////////////////////////////////////////
                }

                c_simd.store(&this->c_[b_offset], alignment);
            };

            for (SizeType b = 0; b < n_block_rows; ++b) {
                f(b);
            }

            if (this->symmetric_) {
                static_assert(std::is_signed<SizeType>::value, "needs to be a signed integer");

                for (SizeType b = n_block_rows - 1; b >= 0; --b) {
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
#endif  // UTOPIA_VC_PROJECTED_BLOCK_GAUSS_SEIDEL_SWEEP_TRANSPOSED_HPP
