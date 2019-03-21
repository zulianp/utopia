#ifndef UTOPIA_PETSC_QUIRKS_HPP
#define UTOPIA_PETSC_QUIRKS_HPP

#include "petscvec.h"    
#include "utopia_Config.hpp"
	

	//Weird name changing in minor versions
	inline PetscErrorCode UtopiaVecScatterCreate(Vec xin, IS ix, Vec yin, IS iy, VecScatter *newctx)
	{
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3,10,4)
		return VecScatterCreate(xin, ix, yin, iy, newctx);
#else
	// #if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3,11,0) || (UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3,10,0) && PETSC_VERSION_RELEASE == 0)
	#if (UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3,10,0) && PETSC_VERSION_RELEASE == 0)
		return VecScatterCreateWithData(xin, ix, yin, iy, newctx);
	#else
		return VecScatterCreate(xin, ix, yin, iy, newctx);
	#endif
#endif
	}

#endif //UTOPIA_PETSC_QUIRKS_HPP
