#ifndef UTOPIA_TPETRA_ROW_VIEW_HPP
#define UTOPIA_TPETRA_ROW_VIEW_HPP
#include <Teuchos_ArrayViewDecl.hpp>

#include "utopia_RowView.hpp"
#include "utopia_Traits.hpp"

#include <cassert>

namespace utopia {

    template <class Tensor, int FILL_TYPE>
    class RowView<Tensor, 2, FILL_TYPE, utopia::TRILINOS> {
    public:
        using SizeType = typename Traits<Tensor>::SizeType;
        using Scalar = typename Traits<Tensor>::Scalar;
        using LocalSizeType = typename Traits<Tensor>::LocalSizeType;

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
        using LocalIndsViewType = typename Traits<Tensor>::Matrix::CrsMatrixType::local_inds_host_view_type;
        using HostViewType = typename Traits<Tensor>::Matrix::CrsMatrixType::values_host_view_type;
#else
        using LocalIndsViewType = Teuchos::ArrayView<const LocalSizeType>;
        using HostViewType = Teuchos::ArrayView<const Scalar>;
#endif

        using RCPMapType = typename Tensor::RCPMapType;

        inline RowView(const Tensor &t, const SizeType row) : t_(t), offset_(0) {
            auto impl = raw_type(t_);
            assert(impl->isLocallyIndexed());
            auto rr = row_range(t);
            impl->getLocalRowView(row - rr.begin(), cols_, values_);
            offset_ = impl->getDomainMap()->getMinGlobalIndex();
            col_map_ = impl->getColMap();
        }

        inline ~RowView() {}

        inline std::size_t n_values() const { return cols_.size(); }

        inline SizeType col(const int index) const {
            assert(static_cast<std::size_t>(index) < n_values());

            auto ret = col_map_->getGlobalElement(cols_[index]);
            assert(ret < SizeType(size(t_).get(1)));
            return ret;
        }

        inline Scalar get(const int index) const {
            assert(static_cast<std::size_t>(index) < n_values());
            return values_[index];
        }

    private:
        const Tensor &t_;
        SizeType offset_;
        LocalIndsViewType cols_;
        HostViewType values_;
        RCPMapType col_map_;
    };
}  // namespace utopia

#endif  // UTOPIA_TPETRA_ROW_VIEW_HPP
