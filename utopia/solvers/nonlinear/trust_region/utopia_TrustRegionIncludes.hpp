#ifndef UTOPIA_TR_INCLUDES_HPP
#define UTOPIA_TR_INCLUDES_HPP


#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_TR_base.hpp"
#include "utopia_TR_NormalEquation.hpp" 

#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_ActiveSetTRSubproblem.hpp"
#include "utopia_TrustRegion.hpp"




#include "utopia_TR_Factory.hpp"


#ifdef WITH_PETSC
	#include "utopia_PETScKSP_TR.hpp"
#endif

#endif //UTOPIA_TR_INCLUDES_HPP