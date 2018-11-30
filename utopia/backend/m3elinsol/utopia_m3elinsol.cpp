#include "utopia_m3elinsol_impl.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos.hpp"
#endif //WITH_TRILINOS

#ifdef WITH_PETSC
#include "utopia_petsc.hpp"
#endif //WITH_PETSC

#ifdef WITH_BLAS
#include "utopia_blas.hpp"
#endif //WITH_BLAS

//explicit instantiations

namespace utopia {

#ifdef WITH_PETSC
	template class ASPAMG<DSMatrixd, DVectord>;
#endif //WITH_PETSC

#ifdef WITH_BLAS

	template class ASPAMG<CRSMatrixd, Vectord>;
#endif //WITH_BLAS

#ifdef WITH_TRILINOS

	template class ASPAMG<TSMatrixd, TVectord>;
#endif //WITH_TRILINOS

}
