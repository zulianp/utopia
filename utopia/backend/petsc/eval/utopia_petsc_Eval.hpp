#ifndef UTOPIA_EVAL_PETSC_HPP
#define UTOPIA_EVAL_PETSC_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Eval_Factory.hpp"
#include "utopia_petsc_Eval_Inverse.hpp"
#include "utopia_petsc_Eval_KroneckerProduct.hpp"
#include "utopia_petsc_Eval_Parallel.hpp"
#include "utopia_petsc_Eval_Rename.hpp"
// #include "utopia_petsc_Eval_Chop.hpp" //UNCOMMENT ME TO USE ALENA'S VERSION
#include "utopia_petsc_EvalDotVecVecs.hpp"
#include "utopia_petsc_EvalMatGetCol.hpp"
#include "utopia_petsc_Eval_Blocks.hpp"
#include "utopia_petsc_Eval_DotOpDot.hpp"
#include "utopia_petsc_Eval_Misc.hpp"
#include "utopia_petsc_Eval_Monitor.hpp"
#include "utopia_petsc_Eval_Multiply.hpp"
#include "utopia_petsc_Eval_NZZXRow.hpp"
#include "utopia_petsc_Eval_Residual.hpp"
#include "utopia_petsc_Eval_VecUniqueSortSerial.hpp"

#ifdef WITH_SLEPC
#include "utopia_petsc_Eval_Cond.hpp"
#endif  // WITH_SLEPC

/*! @file
 * Petsc language extensions
 */

namespace utopia {}

#endif  // UTOPIA_EVAL_PETSC_HPP
