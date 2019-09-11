#ifndef UTOPIA_BLAS_ROW_VIEW_HPP
#define UTOPIA_BLAS_ROW_VIEW_HPP

#include "utopia_Traits.hpp"
#include "utopia_blas_Types.hpp"

namespace utopia {

    template<class Tensor>
    class RowView<Tensor, 2, FillType::DENSE, utopia::BLAS> {
    public:
        typedef typename utopia::Traits<Tensor>::SizeType SizeType;
        typedef typename utopia::Traits<Tensor>::Scalar Scalar;

        inline RowView(Tensor &t, const SizeType row)
        : t_(t), row_(row)
        {}

        inline ~RowView() {}

        inline SizeType n_values() const
        {
            return t_.implementation().cols();
        }

        inline SizeType col(const SizeType index) const
        {
            return index;
        }

        inline Scalar get(const SizeType index) const
        {
            return t_.get(row_, index);
        }

    private:
        Tensor &t_;
        SizeType row_;
    };

    // template<>
    // class RowView<CRSMatrixd, 2, FillType::SPARSE, utopia::BLAS> {
    // public:
    //     typedef typename utopia::Traits<CRSMatrixd>::SizeType SizeType;
    //     typedef typename utopia::Traits<CRSMatrixd>::Scalar Scalar;

    //     inline RowView(CRSMatrixd &t, const SizeType row)
    //     : t_(t),
    //       row_(row),
    //       n_cols_(t.implementation().rowptr()[row+1] - t.implementation().rowptr()[row])
    //     {}

    //     inline ~RowView() {}

    //     inline SizeType n_values() const
    //     {
    //         return n_cols_;
    //     }

    //     inline SizeType col(const SizeType index) const
    //     {
    //         assert(index < n_values());
    //         return t_.implementation().colindex()[t_.implementation().rowptr()[row_] + index];
    //     }

    //     inline Scalar get(const SizeType index) const
    //     {
    //         assert(index < n_values());
    //         return t_.implementation().entries()[t_.implementation().rowptr()[row_] + index];
    //     }

    // private:
    //     CRSMatrixd &t_;
    //     SizeType row_;
    //     SizeType n_cols_;
    // };


    // template<>
    // class RowView<const CRSMatrixd, 2, FillType::SPARSE, utopia::BLAS> {
    // public:
    //     typedef typename utopia::Traits<CRSMatrixd>::SizeType SizeType;
    //     typedef typename utopia::Traits<CRSMatrixd>::Scalar Scalar;

    //     inline RowView(const CRSMatrixd &t, const SizeType row)
    //     : t_(t),
    //       row_(row),
    //       n_cols_(t.implementation().rowptr()[row + 1] - t.implementation().rowptr()[row])
    //     {}

    //     inline ~RowView() {}

    //     inline SizeType n_values() const
    //     {
    //         return n_cols_;
    //     }

    //     inline SizeType col(const SizeType index) const
    //     {
    //         assert(index < n_values());
    //         return t_.implementation().colindex()[t_.implementation().rowptr()[row_] + index];
    //     }

    //     inline Scalar get(const SizeType index) const
    //     {
    //         assert(index < n_values());
    //         return t_.implementation().entries()[t_.implementation().rowptr()[row_] + index];
    //     }

    // private:
    //     const CRSMatrixd &t_;
    //     SizeType row_;
    //     SizeType n_cols_;
    // };

}

#endif //UTOPIA_BLAS_ROW_VIEW_HPP
