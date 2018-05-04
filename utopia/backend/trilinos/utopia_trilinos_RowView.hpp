#ifndef UTOPIA_TPETRA_ROW_VIEW_HPP
#define UTOPIA_TPETRA_ROW_VIEW_HPP
#include <Teuchos_ArrayViewDecl.hpp>

namespace utopia {

	template<class Tensor, int FILL_TYPE>
	class RowView<Tensor, 2, FILL_TYPE, utopia::TRILINOS> {
	public:
		typedef typename Tensor::Implementation::global_ordinal_type global_ordinal_type;
		typedef typename Tensor::Implementation::Scalar Scalar;

		inline RowView(const Tensor &t, const global_ordinal_type row)
		: t_(t), offset_(0)
		{
			if(t_.implementation().implementation().isGloballyIndexed()) {
				t_.implementation().implementation().getGlobalRowView(row, cols_, values_);
			} else {
				assert(t_.implementation().implementation().isLocallyIndexed());
				t_.implementation().implementation().getLocalRowView(row, cols_, values_);
				offset_ = t_.implementation().implementation().getColMap()->getMinGlobalIndex();
			}
		}

		inline ~RowView()
		{}

		inline std::size_t n_values() const
		{
			return cols_.size();
		}

		inline global_ordinal_type col(const int index) const
		{
			assert(index < n_values());
			return cols_[index] + offset_;
		}

		inline Scalar get(const int index) const
		{
			assert(index < n_values());
			return values_[index];
		}

	private:
		const Tensor &t_;
		global_ordinal_type offset_;
		Teuchos::ArrayView<const global_ordinal_type> cols_;
		Teuchos::ArrayView<const Scalar> values_;
	};
}

#endif //UTOPIA_TPETRA_ROW_VIEW_HPP
