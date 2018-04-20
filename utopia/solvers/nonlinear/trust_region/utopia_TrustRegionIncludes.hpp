#ifndef UTOPIA_TR_INCLUDES_HPP
#define UTOPIA_TR_INCLUDES_HPP

// subproblems 
#include "utopia_TRSubproblem.hpp"
#include "utopia_CauchyPoint.hpp"
#include "utopia_Dogleg.hpp"
#include "utopia_SteihaugToint.hpp"
#include "utopia_ActiveSetTRSubproblem.hpp"
// there is some problem with includes inside of petsc 
// base classes 
#include "utopia_TRBase.hpp"
#include "utopia_TRBoxBase.hpp"

// tr-based minimization algorithms 
#include "utopia_TrustRegion.hpp"
#include "utopia_TRNormalEquation.hpp" 
#include "utopia_TrustRegionVariableBound.hpp"


// factory 
#include "utopia_TrustRegionFactory.hpp"


#endif //UTOPIA_TR_INCLUDES_HPP