#ifndef UTOPIA_TRILINOS_SOLVERS_HPP
#define UTOPIA_TRILINOS_SOLVERS_HPP

#include "utopia_PreconditionedSolver.hpp"
//#include "utopia_trilinos_LinearSolverFactory.hpp"

#ifdef WITH_TRILINOS_BELOS
#include "utopia_Belos_solver.hpp"
#endif

#include "utopia_Amesos2_solver.hpp"

#include "utopia_trilinos_ConjugateGradient.hpp"

#endif  // UTOPIA_TRILINOS_SOLVERS_HPP
