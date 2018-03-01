#ifndef UTOPIA_PETSC_HPP
#define UTOPIA_PETSC_HPP 

#include "utopia_petsc_Backend.hpp"
#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_SerialSparseMatrix.hpp"
#include "utopia_petsc_SerialVector.hpp"
#include "utopia_petsc_SparseMatrix.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Eval.hpp"
#include "utopia_petsc_Eval_Inverse.hpp"
#include "utopia_petsc_Eval_Factory.hpp"
#include "utopia_petsc_RowView.hpp"

#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_TrustRegionFactory.hpp"
#include "utopia_petsc_solvers.hpp"

#include "utopia_petsc_KSPTR.hpp"

namespace utopia {
	void optimize_nnz(DSMatrixd &A);
}

#endif //UTOPIA_PETSC_HPP
