#include "utopia_petsc_KSPWrapper.hpp"
#include "utopia_petsc_Types.hpp"



namespace utopia {
	KSPTypes::KSPTypes():
	ksp_({"bcgs", "cg", "groppcg", "pipecg", "pipecgrr", "fcg", "pipefcg", "gmres", "pipefgmres",   "fgmres",   "lgmres",   "dgmres",   "pgmres", "tcqmr", "ibcgs",   "fbcgs",   "fbcgsr",   "bcgsl", "cgs", "tfqmr", "cr", "pipecr", "lsqr", "preonly", "qcg", "bicg", "minres", "symmlq", "lcd", "python", "gcr", "pipegcr", "tsirm", "cgls"}),
	pc_({"jacobi","sor","lu","bjacobi","eisenstat","ilu","icc","asm","gasm","ksp","cholesky","pbjacobi","mat","hypre", "cp","bfbt","lsc","python","pfmg","syspfmg","redistribute","svd","gamg","bicgstabcusp","ainvcusp","bddc"}), 
	package_({" "})
	{}

	template class KSPWrapper<DSMatrixd, DVectord>;
}
