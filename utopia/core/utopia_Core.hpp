#ifndef UTOPIA_CORE_HPP
#define UTOPIA_CORE_HPP


// to get rid of unused variables in Realease build 
// #define ASSERT_DEBUG

#ifdef ASSERT_DEBUG
#define ASSERT(x) do { (void)sizeof(x);} while (0)
#else
#include <assert.h>
#define ASSERT(x) assert(x)
#endif


#include "utopia_Expressions.hpp"
#include "utopia_AuxiliaryExpressions.hpp"
#include "utopia_ExpressionsParallel.hpp"



/** @defgroup base_functions Base Functions
 *  @brief      Base functions used in utopia programms.
 */


#include "utopia_Backend.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Evaluator.hpp"
#include "utopia_Factory.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Instance.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Each.hpp"
#include "utopia_RowView.hpp"
#include "utopia_Unfold.hpp"
#include "utopia_TreeProperties.hpp"
#include "utopia_ScalarCast.hpp"
#include "utopia_TreeNavigator.hpp"
#include "utopia_Variable.hpp"
#include "utopia_ExprInliner.hpp"
#include "utopia_Eval.hpp"

#include "utopia_Temp.hpp"

#endif //UTOPIA_CORE_HPP
