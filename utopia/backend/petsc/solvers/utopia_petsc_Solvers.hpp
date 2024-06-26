#ifndef UTOPIA_PETSC_SOLVERS_HPP
#define UTOPIA_PETSC_SOLVERS_HPP

#include "utopia_petsc_LinearSolvers.hpp"
#include "utopia_petsc_MultilevelSolvers.hpp"
#include "utopia_petsc_NonlinearSolvers.hpp"
#include "utopia_petsc_Smoothers.hpp"

#ifdef UTOPIA_ENABLE_SLEPC
#include "utopia_petsc_Slepc.hpp"
#endif

#endif  // UTOPIA_PETSC_SOLVERS_HPP