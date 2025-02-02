#ifndef UTOPIA_HPP
#define UTOPIA_HPP

/*! @file utopia.hpp
 *  Include this file for using utopia.
 */

/*! @mainpage
* \tableofcontents
*
*
* @section Utopia a C++ library for large scale scientific computing
* The main idea behind Utopia's design is the separation of the application
layer and the technology used to realize it.
* Following this idea we describe applications through Utopia's EDSL which by
itself does not compute anything.
* In fact, the EDSL exclusively generates a statically typed expression tree.
This expression tree can be used or transformed for different purposes.
* Its most obvious usage is evaluation of the represented expression (e.g.,
@code b = A*x; @endcode), however it can be transformed to fit other purposes
such as expression optimization (e.g.,  @em (A * B) * x, where @em A and @em B
are matrices and @em x a vector, becomes A * ( B * x)), simplification or matrix
differentiation. An example of matrix differentiation is the computation of the
derivative of
@code transpose(x) * A * x + transpose(x) * b @endcode
with respect to @em x which is transformed without any computation to
@code A * x + b @endcode
with compile time decisions.
* The manipulation of the expression tree is still an application layer activity
and no actual computation is yet performed.

</br>
Visit our <a href="https://bitbucket.org/zulianp/utopia/wiki/Home">Wiki</a> for
more information and tutorials.
*/

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Solvers.hpp"
#include "utopia_Utils.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS
// Macros defined in a petsc header are conflicting with Enums defined in
// Trilinos Trilinos inclusions should always preceed petsc inclusions
#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"
#endif  // UTOPIA_ENABLE_TRILINOS

#ifdef UTOPIA_ENABLE_BLAS
#include "utopia_blas.hpp"
#endif  // UTOPIA_ENABLE_BLAS

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"
#endif  // UTOPIA_ENABLE_PETSC

#ifdef UTOPIA_ENABLE_VC
#include "utopia_Vc.hpp"
#endif  // UTOPIA_ENABLE_VC

#include "utopia_polymorphic_LinearSolver.hpp"

#endif  // UTOPIA_HPP
