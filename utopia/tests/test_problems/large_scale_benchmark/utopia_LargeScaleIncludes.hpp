#ifndef UTOPIA_LARGE_SCALE_TEST_FUNCTIONS_HPP
#define UTOPIA_LARGE_SCALE_TEST_FUNCTIONS_HPP

#include "utopia_Bratu1D.hpp"
#include "utopia_MOREBV1D.hpp"
#include "utopia_Poisson1D.hpp"

#include "utopia_Bratu2D.hpp"
#include "utopia_Bratu3D.hpp"
#include "utopia_Membrane2D.hpp"
#include "utopia_NonEllipse2D.hpp"
#include "utopia_Poisson2D.hpp"
#include "utopia_Poisson3D.hpp"

#include "utopia_MultilevelTestProblem1D.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_PetscMultilevelTestProblem.hpp"
#endif

#endif  // UTOPIA_LARGE_SCALE_TEST_FUNCTIONS_HPP
