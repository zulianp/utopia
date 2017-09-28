#ifndef UTOPIA_PETSC_HPP
#define UTOPIA_PETSC_HPP 

#include "utopia_PETScBackend.hpp"
#include "utopia_PETScError.hpp"
#include "utopia_PETScForwardDeclaration.hpp"
#include "utopia_PETScMatrix.hpp"
#include "utopia_PETScSerialSparseMatrix.hpp"
#include "utopia_PETScSerialVector.hpp"
#include "utopia_PETScSparseMatrix.hpp"
#include "utopia_PETScTraits.hpp"
#include "utopia_PETScTypes.hpp"
#include "utopia_PETScVector.hpp"

#include "utopia_Eval_PETSc.hpp"
#include "utopia_Eval_Inverse_PETSc.hpp"
#include "utopia_petsc_RowView.hpp"

#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_TrustRegionFactory.hpp"
#include "utopia_petsc_solvers.hpp"

#include "utopia_PETScLinearGS.hpp"
#include "utopia_PETScCGSmoother.hpp"
#include "utopia_PETScNonLinearGS.hpp"
#include "utopia_PETScFunction.hpp"

#include "utopia_PETScKSPTR.hpp"

#endif //UTOPIA_PETSC_HPP
