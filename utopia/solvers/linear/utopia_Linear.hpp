/*
* @Author: alenakopanicakova
* @Date:   2016-05-10
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-03
*/
#ifndef SOLVER_SOLVER_LINEAR_INCLUDES_HPP
#define SOLVER_SOLVER_LINEAR_INCLUDES_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_DirectSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_ConjugateGradient.hpp"

#ifdef WITH_PETSC
#include "utopia_PETScFactorizations.hpp"
#include "utopia_PETScKSPSolvers.hpp"
#endif //WITH_PETSC

#ifdef WITH_LAPACK
#include "utopia_Lapack.hpp"
#endif //WITH_LAPACK


#ifdef WITH_UMFPACK
#include "utopia_UmfpackLU.hpp"
#endif //WITH_UMFPACK

#include "utopia_LinearSolverFactory.hpp"

#endif