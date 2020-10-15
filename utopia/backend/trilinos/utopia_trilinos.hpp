#ifndef UTOPIA_TRILINOS_HPP
#define UTOPIA_TRILINOS_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#include "utopia_Kokkos_ParallelFor.hpp"
#include "utopia_trilinos_DeviceView.hpp"
#include "utopia_trilinos_Eval_Factory.hpp"
#include "utopia_trilinos_ForwardDeclarations.hpp"
#include "utopia_trilinos_RowView.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_trilinos_Types.hpp"
#include "utopia_trilinos_solvers.hpp"

// optimizations
#include "utopia_trilinos_Eval_RAP.hpp"

// FIXME this should not be necessary
#include "utopia_kokkos_Eval_Binary.hpp"
#include "utopia_kokkos_Eval_MultiReduce.hpp"
#include "utopia_kokkos_Eval_Reduce.hpp"
#include "utopia_kokkos_Eval_Unary.hpp"

#include "utopia_Tpetra_Matrix_impl.hpp"
#include "utopia_Tpetra_Vector_impl.hpp"
#include "utopia_trilinos_DiffController.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"
#include "utopia_trilinos_MaxRowNNZ.hpp"
// FIXME re-introduce later
// #include "utopia_trilinos_Eval_Distance.hpp"

#endif  // UTOPIA_WITH_TRILINOS
#endif  // UTOPIA_TRILINOS_HPP
