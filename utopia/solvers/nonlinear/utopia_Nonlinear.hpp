#ifndef UTOPIA_NON_LINEAR_HPP
#define UTOPIA_NON_LINEAR_HPP 


#include "utopia_Function.hpp"
#include "utopia_GLFunction.hpp"
#include "utopia_FunctionNormalEq.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_Newton.hpp"


#include "utopia_InexactNewton.hpp"

#include "utopia_ConstrainedIncludes.hpp"
#include "utopia_TrustRegionIncludes.hpp"
#include "utopia_LineSearchIncludes.hpp"


#include "utopia_NonlinearSolverFactory.hpp"


#ifdef WITH_PETSC
	#include "utopia_PETScFunction.hpp"
#endif //WITH_PETSC	

#endif //UTOPIA_NON_LINEAR_HPP

