#ifndef UTOPIA_TPETRA_ROW_VIEW_HPP
#define UTOPIA_TPETRA_ROW_VIEW_HPP
#include <Teuchos_ArrayViewDecl.hpp>

#include "utopia_Traits.hpp"
#include <cassert>

namespace utopia {

	template<class Tensor, int FILL_TYPE>
	class RowView<Tensor, 2, FILL_TYPE, utopia::TRILINOS> {
	public:
		typedef typename Tensor::Implementation::GO GO;
		typedef typename Tensor::Implementation::Scalar Scalar;

		inline RowView(const Tensor &t, const GO row, const bool force_local_view = true)
		: t_(t), offset_(0)
		{
			assert(t_.implementation().implementation().isLocallyIndexed());
			auto rr = row_range(t);
			t_.implementation().implementation().getLocalRowView(row - rr.begin(), cols_, values_);
			offset_ = t_.implementation().implementation().getDomainMap()->getMinGlobalIndex();
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

			auto ret = t_.implementation().implementation().getColMap()->getGlobalElement(cols_[index]);
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
	};
}

#endif //UTOPIA_TPETRA_ROW_VIEW_HPP
