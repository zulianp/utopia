#ifndef UTOPIA_PETSC_EACH_HPP
#define UTOPIA_PETSC_EACH_HPP

#include "utopia_petsc_Types.hpp"
#include "utopia_petsc_RowView.hpp"

namespace utopia {

	template<class Fun>
	inline void each_apply(DSMatrixd &mat, Fun fun) 
	{	
		//FIXME Very innefficient but bust find out other way
		DSMatrixd mat_copy = mat;

		Write<DSMatrixd> w(mat);
		each_read(mat_copy, [&mat, &fun](const PetscInt i, const PetscInt j, const PetscScalar value) {
			mat.set(i, j, fun(value));
		});
	}
}

#endif //UTOPIA_PETSC_EACH_HPP
