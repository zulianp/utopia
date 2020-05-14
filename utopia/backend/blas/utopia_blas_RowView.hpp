#ifndef UTOPIA_BLAS_ROW_VIEW_HPP
#define UTOPIA_BLAS_ROW_VIEW_HPP

#include "utopia_Traits.hpp"
#include "utopia_blas_Types.hpp"

namespace utopia {

    template<class Tensor>
    class RowView<Tensor, 2, FillType::DENSE, utopia::BLAS> {
    public:
        using SizeType = typename utopia::Traits<Tensor>::SizeType;
        using Scalar = typename utopia::Traits<Tensor>::Scalar;

        inline RowView(Tensor &t, const SizeType row)
        : t_(t), row_(row)
        {}

        inline ~RowView() = default;

        inline SizeType n_values() const
        {
            return t_.cols();
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

}

#endif //UTOPIA_BLAS_ROW_VIEW_HPP
