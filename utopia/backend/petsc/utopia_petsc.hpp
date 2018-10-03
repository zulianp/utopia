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
#include "utopia_petsc_Eval_Residual.hpp"
#include "utopia_petsc_Eval_DotDivDot.hpp"
#include "utopia_petsc_Eval_Blocks.hpp"
#include "utopia_petsc_RowView.hpp"
#include "utopia_petsc_EvalDotVecVecs.hpp"

#include "utopia_petsc_LinearSolverFactory.hpp"
#include "utopia_petsc_TrustRegionFactory.hpp"
#include "utopia_petsc_solvers.hpp"

#include "utopia_petsc_Newton.hpp"
#include "utopia_petsc_TaoTRQP.hpp"



#ifdef WITH_SLEPC
	#include "utopia_petsc_Slepc.hpp"
	#include "utopia_petsc_MoreSorensenEigen.hpp"
#endif


// very much experimental files for the moment 
#include "utopia_petsc_SNES.hpp"
#include "utopia_petsc_build_ksp.hpp"


namespace utopia {
	void optimize_nnz(DSMatrixd &A);
	bool is_diagonally_dominant(const DSMatrixd &A);
}

#endif //UTOPIA_PETSC_HPP
