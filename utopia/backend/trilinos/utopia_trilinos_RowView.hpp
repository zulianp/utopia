#ifndef UTOPIA_TPETRA_ROW_VIEW_HPP
#define UTOPIA_TPETRA_ROW_VIEW_HPP
#include <Teuchos_ArrayViewDecl.hpp>

#include "utopia_Traits.hpp"
#include "utopia_RowView.hpp"

#include <cassert>

namespace utopia {

	template<class Tensor, int FILL_TYPE>
	class RowView<Tensor, 2, FILL_TYPE, utopia::TRILINOS> {
	public:
		typedef typename Tensor::Implementation::GO GO;
		typedef typename Tensor::Implementation::Scalar Scalar;
		using MapPtrT = typename Tensor::Implementation::rcp_map_type;

		inline RowView(const Tensor &t, const GO row, const bool force_local_view = true)
		: t_(t), offset_(0)
		{
			auto impl = raw_type(t_);
			assert(impl->isLocallyIndexed());
			auto rr = row_range(t);
			impl->getLocalRowView(row - rr.begin(), cols_, values_);
			offset_  = impl->getDomainMap()->getMinGlobalIndex();
			col_map_ = impl->getColMap();
		}

		inline ~RowView()
		{}

		inline std::size_t n_values() const
		{
			return cols_.size();
		}

		inline GO col(const int index) const
		{
			assert(index < n_values());

			auto ret = col_map_->getGlobalElement(cols_[index]);
			assert(ret < size(t_).get(1));
			return ret;
		}

		inline Scalar get(const int index) const
		{
			assert(index < n_values());
			return values_[index];
		}

	private:
		const Tensor &t_;
		GO offset_;
		Teuchos::ArrayView<const GO> cols_;
		Teuchos::ArrayView<const Scalar> values_;
		MapPtrT col_map_;
	};
}

#endif //UTOPIA_TPETRA_ROW_VIEW_HPP
