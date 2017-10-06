#ifndef UTOPIA_CONSTRAINTED_HPP
#define UTOPIA_CONSTRAINTED_HPP 

#include "utopia_Base.hpp"
#include "utopia_NonlinSemismoothNewton.hpp"
#include "utopia_MultigridConstrained.hpp"

#include "utopia_FunctionBoxConstrained.hpp"

#include "utopia_MPRGP.hpp"
#include "utopia_SemismoothNewton.hpp"

#ifdef WITH_PETSC
#include "utopia_PETScSemismoothNewton.hpp"	
#endif

#endif //UTOPIA_CONSTRAINTED_HPP

