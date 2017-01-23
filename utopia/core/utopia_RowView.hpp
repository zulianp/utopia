#ifndef UTOPIA_ROW_VIEW_HPP
#define UTOPIA_ROW_VIEW_HPP 

#include "utopia_Base.hpp"

namespace utopia {

	template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Tensor::FILL_TYPE, int BackendType = Traits<Tensor>::Backend>
	class RowView {};


#ifdef WITH_PETSC
	template<class Tensor, int FILL_TYPE>
	class RowView<Tensor, 2, FILL_TYPE, utopia::PETSC> {
	public:
		//FIXME move petsc specific stuff to backend
		
		inline RowView(Tensor &t, const PetscInt row)
		: t_(t), row_(row)
		{
			MatGetRow(raw_type(t_), row_, &n_values_, &cols_, &values_);
		}

		inline ~RowView()
		{
			MatRestoreRow(raw_type(t_), row_, &n_values_, &cols_, &values_);
		}
		
		inline PetscInt n_values() const
		{
			return n_values_;
		}

		inline PetscInt get_col_at(const PetscInt index) const
		{
			assert(index < n_values_);
			return cols_[index];
		}

		inline PetscScalar get_value_at(const PetscInt index) const
		{
			assert(index < n_values_);
			return values_[index];
		}

		// inline void set_value_at(const PetscInt index, const PetscScalar value) const
		// {
		// 	assert(index < n_values_);
		// 	values_[index] = value;
		// }

	private:
		Tensor &t_;
		PetscInt row_;
		PetscInt n_values_;
		const PetscInt 	  * cols_;
		const PetscScalar * values_;
	};
#endif //WITH_PETSC	
}


#endif //UTOPIA_ROW_VIEW_HPP

