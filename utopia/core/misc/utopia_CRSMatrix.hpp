#ifndef UTOPIA_CRS_MATRIX_HPP
#define UTOPIA_CRS_MATRIX_HPP

#include "utopia_Base.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

#include <vector>

namespace utopia {

    template <typename T>
    class Traits<std::vector<T> > {
    public:
        using ValueType = T;
    };

    template <class ScalarView, class IndexView, int BlockSize_ = 1>
    class CRSMatrix {
    public:
        using Scalar = typename utopia::Traits<ScalarView>::ValueType;
        using SizeType = typename utopia::Traits<IndexView>::ValueType;
        static const int BlockSize = BlockSize_;
        static const int BlockSize2 = BlockSize * BlockSize;

        class SparseRowView {
        public:
            UTOPIA_INLINE_FUNCTION Scalar &value(const SizeType i) { return values_[i]; }

            UTOPIA_INLINE_FUNCTION const Scalar &value(const SizeType i) const { return values_[i]; }

            UTOPIA_INLINE_FUNCTION Scalar *block(const SizeType i) const { return &values_[i * BlockSize2]; }
            UTOPIA_INLINE_FUNCTION SizeType n_blocks() const { return length / BlockSize2; }

            UTOPIA_INLINE_FUNCTION const SizeType &colidx(const SizeType i) const { return colidx_[i]; }

            UTOPIA_INLINE_FUNCTION SizeType &colidx(const SizeType i) { return colidx_[i]; }

            UTOPIA_INLINE_FUNCTION SparseRowView(Scalar *values, SizeType *colidx, const SizeType &length)
                : length(length), values_(values), colidx_(colidx) {}

            const SizeType length;

        private:
            Scalar *values_;
            SizeType *colidx_;
        };

        UTOPIA_INLINE_FUNCTION CRSMatrix(IndexView row_ptr, IndexView colidx, ScalarView values, const SizeType n_cols)
            : row_ptr_(row_ptr), colidx_(colidx), values_(values), n_cols_(n_cols) {}

        UTOPIA_INLINE_FUNCTION SparseRowView row(const SizeType i) {
            return SparseRowView(&values_[row_ptr_[i] * BlockSize2],
                                 &colidx_[row_ptr_[i]],
                                 (row_ptr_[i + 1] - row_ptr_[i]) * BlockSize2);
        }

        UTOPIA_INLINE_FUNCTION CRSMatrix() = default;

        UTOPIA_INLINE_FUNCTION SizeType nnz() const { return values_.size(); }
        UTOPIA_INLINE_FUNCTION SizeType cols() const { return n_cols_; }
        UTOPIA_INLINE_FUNCTION SizeType rows() const { return row_ptr_.size() - 1; }

        UTOPIA_INLINE_FUNCTION IndexView &row_ptr() { return row_ptr_; }
        UTOPIA_INLINE_FUNCTION IndexView &colidx() { return colidx_; }
        UTOPIA_INLINE_FUNCTION ScalarView &values() { return values_; }

        UTOPIA_INLINE_FUNCTION const IndexView &row_ptr() const { return row_ptr_; }
        UTOPIA_INLINE_FUNCTION const IndexView &colidx() const { return colidx_; }
        UTOPIA_INLINE_FUNCTION const ScalarView &values() const { return values_; }

        UTOPIA_INLINE_FUNCTION const Scalar *block(const SizeType i) const { return &values_[i * BlockSize2]; }
        UTOPIA_INLINE_FUNCTION Scalar *block(const SizeType i) { return &values_[i * BlockSize2]; }

        UTOPIA_INLINE_FUNCTION void set_cols(const SizeType n_cols) { n_cols_ = n_cols; }

    private:
        IndexView row_ptr_;
        IndexView colidx_;
        ScalarView values_;
        SizeType n_cols_{-1};
    };

    template <typename S, typename I, int BlockSize>
    class Traits<CRSMatrix<S, I, BlockSize> > {
    public:
        static const int Backend = HOMEMADE;
        using Scalar = typename Traits<I>::ValueType;
        using SizeType = typename Traits<I>::ValueType;
    };

}  // namespace utopia

#endif  // UTOPIA_CRS_MATRIX_HPP
