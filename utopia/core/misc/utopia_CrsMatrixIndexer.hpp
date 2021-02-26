#ifndef UTOPIA_CRS_MATRIX_INDEXER_HPP
#define UTOPIA_CRS_MATRIX_INDEXER_HPP

#include <algorithm>
#include <vector>

namespace utopia {

    template <typename SizeType>
    class CrsDiagIndexer {
    public:
        void init(SizeType n, const SizeType *row_ptr, const SizeType *colidx) {
            idx_.resize(n);
            this->row_ptr_ = row_ptr;
            this->colidx_ = colidx;

            for (SizeType r = 0; r < n; ++r) {
                const SizeType row_end = row_ptr_[r + 1];

                SizeType k;
                bool has_diag = false;
                for (k = row_ptr_[r]; k < row_end; ++k) {
                    if (colidx_[k] == r) {
                        has_diag = true;
                        break;
                    }
                }

                assert(has_diag);
                UTOPIA_UNUSED(has_diag);

                idx_[r] = k;
            }
        }

        SizeType find(const SizeType i, const SizeType j) const {
            if (i == j) return idx_[i];

            if (j < i) {
                return find_before_diag(i, j);
            } else {
                return find_after_diag(i, j);
            }
        }

        inline SizeType diag(const SizeType i) const { return idx_[i]; }

    private:
        std::vector<SizeType> idx_;
        const SizeType *row_ptr_;
        const SizeType *colidx_;

        auto find_before_diag(const SizeType i, const SizeType j) const -> SizeType {
            SizeType search_end = idx_[i];
            SizeType row_begin = row_ptr_[i];

            assert(j <= i);
            assert(row_begin <= search_end);
            assert(std::is_sorted(colidx_ + row_begin, colidx_ + search_end));

            auto s_begin = colidx_ + row_begin;
            auto s_end = colidx_ + search_end;

            const SizeType *it = std::lower_bound(s_begin, s_end, j);
            const SizeType ret = (it == s_end || j != *it) ? SizeType(-1) : SizeType(it - colidx_);

#ifndef NDEBUG
            SizeType actual_ret = -1;
            for (SizeType k = row_begin; k < search_end; ++k) {
                if (colidx_[k] == j) {
                    actual_ret = k;
                    break;
                }
            }

            assert(ret == actual_ret);
#endif  // NDEBUG

            return ret;
        };

        auto find_after_diag(const SizeType i, const SizeType j) const -> SizeType {
            SizeType search_begin = idx_[i];
            SizeType row_end = row_ptr_[i + 1];

            assert(i <= j);
            assert(search_begin <= row_end);
            assert(std::is_sorted(colidx_ + search_begin, colidx_ + row_end));

            auto first = std::lower_bound(colidx_ + search_begin, colidx_ + row_end, j);
            if (first == (colidx_ + row_end) || j != *first) return -1;
            return (first - colidx_);
        };
    };

    template <typename SizeType>
    class CrsTransposeIndexer {
    public:
        void init(const SizeType n, SizeType nnz, const SizeType *row_ptr, const SizeType *colidx) {
            idx_.resize(nnz);

            std::fill(std::begin(idx_), std::end(idx_), -1);

            for (SizeType i = 0; i < n; ++i) {
                const auto row_begin = row_ptr[i];
                const auto row_end = row_ptr[i + 1];

                for (SizeType k = row_begin; k < row_end; ++k) {
                    if (idx_[k] != -1) continue;

                    const auto j = colidx[k];

                    if (i == j) {
                        idx_[k] = k;
                        continue;
                    }

                    const SizeType row_begin_tr = row_ptr[j];
                    const SizeType row_end_tr = row_ptr[j + 1];

                    const SizeType *it = std::lower_bound(colidx + row_begin_tr, colidx + row_end_tr, j);
                    SizeType k_tr = -1;

                    if (it != (colidx + row_end_tr) && j == *it) {
                        k_tr = (it - colidx);
                    }

                    idx_[k] = k_tr;
                    idx_[k_tr] = k;
                }
            }
        }

        inline SizeType transpose(const SizeType k) const {
            assert(k < SizeType(idx_.size()));
            return idx_[k];
        }

    private:
        std::vector<SizeType> idx_;
    };

}  // namespace utopia

#endif  // UTOPIA_CRS_MATRIX_INDEXER_HPP
