#ifndef UTOPIA_HPP
#define UTOPIA_HPP

/*! @file utopia.hpp
 *  Include this file and only this file for using utopia.
 */

/*! @mainpage
* \tableofcontents
*
*
* @section intro A quick introduction to utopia
* The main idea behind Utopia's design is the separation of the application layer and the technology used to realize it.
* Following this idea we describe applications through Utopia's EDSL which by itself does not compute anything.
* In fact, the EDSL exclusively generates a statically typed expression tree. This expression tree can be used or
transformed for different purposes.
* Its most obvious usage is evaluation of the represented expression (e.g., @code b = A*x; @endcode), however it can be
transformed to fit other purposes such as expression optimization (e.g.,  @em (A * B) * x, where @em A and @em B are
matrices and @em x a vector, becomes A * ( B * x)), simplification or matrix differentiation. An example of matrix
differentiation is the computation of the derivative of
@code transpose(x) * A * x + transpose(x) * b @endcode
with respect to @em x which is transformed without any computation to
@code A * x + b @endcode
with compile time decisions.
* The manipulation of the expression tree is still an application layer activity and no actual computation is yet
performed.
*
* With Utopia the evaluation of an expression can be performed in different ways. We identify three main ways: @em
expression to function mapping, @em code @em generation, and @em expression @em templates.
*
* Expression to function mapping
* This type of evaluation is performed by mapping expressions, which are the nodes of the expression tree, to functions
of specific back-ends.
* For instance, the expression @code Vector v = values(n, 0.1); @endcode in our blas based back-end is simply mapped to
a function such as
* @code
* void build(int n, double val, std::vector<double> &v)
* {
* 	v.resize(n);
* 	std::fill(v.begin(), v.end(), val);
* }
* @endcode
* which constructs a vector of length @em n with entries equal to 0.1. The same can be said for all single expressions.
However, more complex composite expression can be pattern-matched and mapped to specific back-end calls; for instance
the expression @em y= alpha * x + y, where @em x, @em y are vectors and @em alpha is a scalar, is mapped to the blas
function @em axpy.
*
* The class in the Utopia library responsible to dispatch expressions or composite expressions to back-end calls is
called @em Evaluator.
*
* In this category of back-ends falls a Petsc based back-end which is currently the most developed in Utopia. However,
Utopia also has a custom made back-end based on the STL, blas, and lapack which is mostly available for learning
purposes.
*
* @section codegen Code generation
* We are currently developing an evaluator/code-generator which generates OpenCL programs, compiles, and runs them
following a just-in-time (JIT) approach (see also ViennaCL). It is still in an experimental state but it aims at very
fast in-lined computations. Statically typed expressions allows generating and compiling specific expression subtrees
only once for each runtime.
*
* @section edsl Utopia EDSL and memory access
* The design of Utopia EDSL is inspired by MatLAB and Eigen and strives for simplicity. However, due to the fact that
Utopia also targets large scale computations, some specific aspects need to be explicitly handled.
* The main aspect is memory location. Although, Utopia provides a certain degree of transparency (as in Petsc) it
requires that things are handled in a memory conscious manner. Independent of which back-end is being used objects and
functions handling ranges and access have to be used always. Ranges allow to iterate over elements of matrices or
vectors which are available in current host address space. For instance, a supercomputer has multiple compute nodes with
dedicated memory (which means that is distributed), hence parts of the data accessible to one node is not directly
accessible to another one. For this example ranges allow to iterate exclusively over the elements that are available
within each node separately.
*
* For a vector @em v its range can be accessed as @em range(v).
* For a matrix @em m its ranges can be accessed as @em row_range(m) and @em col_range(m).
* The utopia::Range class provides the methods @em begin and @em end which are not iterators but just begin and end of
the range.
* See \ref ranges for more details.
*
* In order to read or write from and object we need to acquire its lock and release it when we are done.
* Locks are basic concepts of concurrent programming and are usually used to synchronize access to a shared object.
* We use locks, such as utopia::Read and utopia::Write, in order to provide memory access transparency.
* For instance, the memory region we want to access might not be directly accessible because it is GPU memory hence for
performing CPU operations we need to copy it to the CPU memory.
* When writing in a distributed data-structure such as a distributed matrix the entries we are writing might not be
directly physically accessible, hence a write lock allows communicating all the data at once at the end of the writing
procedure instead of doing each time, which would result in a very slow program.
*
* @code
* Matrix m = zeros(n, n);
* { //<---- Here we create a scope for the lock w.
*
*   //The lock of m is aquired with writing rights,
*   Write<Matrix> w(m);
*   //hence we can add or set any element of m.
*   m.set(0, 0, 1.0);
*
* } //<---- Once we exit this the lock on m is released.
* @endcode
*
* @section devstate Current state of development
* Is currently an experimental library but with very active development.
*
* @section codeexamples Base code examples
*
* In this section we illustrate the basic usage of the utopia EDSL. The examples are not back-end specific hence we use
generic type names for matrices
* (i.e, @em SparseMatrix and @em DenseMatrix) and vectors (i.e, @em Vector)
* which can be directly replaced by specific types
* (e.g., utopia::PetscMatrix, utopia::PetscMatrix and utopia::PetscVector for using the Petsc back-end) or used in a
generic way for algorithms that are not back-end specific.
*
*
* @code
*     // n x n sparse matrix
*     SizeType n = 100;
*     SizeType max_entries_x_row = 3;
*     SparseMatrix m = sparse(n, n, max_entries_x_row);
*
*     {	//Beginning of write lock scope
*         Write<SparseMatrix> w(m);
*         Range r = row_range(m);
*
*         for(SizeType i = r.begin(); i != r.end(); ++i) {
*             if(i > 0) {
*                 m.add(i, i - 1, -1.0);
*             }
*
*             if(i < n-1) {
*                 m.add(i, i + 1, -1.0);
*             }
*
*             m.add(i, i, 2.0);
*         }
*     } 	//End of write lock scope
* @endcode
*
* This example makes use of a factory function called sparse. Utopia provides several factory functions which are
* listed in @ref factory
*
* @section mainclasses Main classes
*	- utopia::Expression is the main abstraction in utopia. In fact, everything is an Expression. More specifically,
*	everything is an abstract (or symbolic) Expression until it is assigned to an utopia::Wrapper. When assigned,
*	the specific backend will evaluate that expression.
* 	- utopia::Wrapper is the main abstraction for tensors. All concrete operations are made on objects of this type.
*	  See the follwing subclasses of utopia::Wrapper for backend specific types:
* 	- petsc wrapper types
*		- utopia::PetscMatrix
*		- utopia::PetscVector
* 	- blas wrapper types
* 		- utopia::BlasMatrix
*		- utopia::BlasVector
*
*
* @section mainfunctions Main functions
* Most typical algebraic expressions and operations are available via operator overloading (* - + / -= *= +=)
* the same way you would do it with Matlab or Eigen. Other operations are available through a functional
* style like interface:
* 	- Permutations
*		- utopia::transpose()
*	- Global factory
*		- utopia::zeros()
*		- utopia::identity()
*		- utopia::values()
*		- utopia::sparse()
*		- utopia::dense()
* 	- Local factory (usefull if your program needs to interface with external software which provides its own data
distribution) *		- utopia::local_identity() *		- utopia::local_zeros() *		-
utopia::local_identity() *		- utopia::local_values() *		- utopia::local_sparse() *	-
Transforms *		- utopia::abs() *		- utopia::sqrt() *		- utopia::pow2() *	-
Reductions *		- utopia::dot() *		- utopia::norm1() *		- utopia::norm2()
*		- utopia::norm_infty()
*		- utopia::sum()
*		- utopia::trace()
*	- Input/output
*		- utopia::read()
*		- utopia::write()
*		- utopia::disp()
*	- Interoperability
*		- utopia::convert()
*		- utopia::raw_type()
*	- Profiling, debugging and user interface
*		- utopia::monitor()
*	- Structural and numerical queries
*		- utopia::size()
*		- utopia::local_size()
*		- utopia::empty()
*		- utopia::range()
*		- utopia::row_range()
*		- utopia::col_range()
*		- utopia::approxeq()
*	- Solving
*		- utopia::solve()
*	- Products
*		- utopia::e_mul()
*		- utopia::outer()
* 	- Utilities
*		- utopia::make_ref()
*
*
* @section solvers Solvers
*	@subsection linear_solvers Linear Solvers
*	@subsubsection precondotioners Precondotioners
*	@subsection nonlinear_solvers Nonlinear Solvers
*	@subsubsection function_description Function
*
*
* @section examples Examples
* @subsection exbasic Basic operations
*
* Setting and getting entries from tensors.
* \snippet tests/utopia_PetscTest.cpp Read write matrix
*
* In place operations (in this case we are using blas types, but it is the same with petsc types).
* \snippet tests/utopia_BlasTest.cpp in place operations (blas)
*
* Example with axpy operation and norms (using petsc types).
* \snippet tests/utopia_PetscTest.cpp axpy (petsc)
*
* Example of input/output operations (using petsc types)
* \snippet tests/utopia_PetscTest.cpp Input and output (petsc)
*
* Views and copy/redistribute a global block of a matrix (using petsc types).
* \snippet tests/utopia_PetscTest.cpp Global views
*
*
* @subsection exsolve Solving systems
* Solving a non-linear system for finding the minimum of a function.
* \snippet tests/utopia_SolverTest.cpp NL solve example
* \snippet tests/utopia_SolverTest.cpp MG solve example
*
*
* @section coding_standards Coding Standard
*  The following coding standards should be followed, when programming Utopia modules and classes:
*  - Class names
* 		- must begin with upcase letters
* 		- use upper case letters as word separator, e.g. BlockJacobi
*  - Method names
*  	    - must begin with lowcase letters
*  	    - use underscore as word separator, e.g. local_zeros
*  - Comment your program
*
*
* @section extending_utopia Extending Utopia
* @subsection extending_utopia_opt Optimizations
* Utopia allows for modular extensions for evaluation optimization which is then applied transparently to the whole
code.
* Consider the following example operation
* @code
* Matrix A = ...;
* Matrix P = ...;
* Matrix m = transpose(P) * A * P;
* @endcode
* This particular operation can be optimized by considering it as a whole for a specific backend, the petsc backend. We
can catch such operations by
* matching the pattern of the expression-tree and by performing some simple checks as follow:
* \snippet core/eval/petsc/utopia_Eval_Petsc.hpp pattern matching and optimizations
* If you want to contribute to utopia or you have some functionality that you want to add (e.g., implemented with pestc
types)
* write us at  <a href="mailto:patrick.zulian@gmail.com"> patrick.zulian@gmail.com </a>.
*
* @section bugs Bug reports
* If you encounter a bug please provide your code snippet responsible for the bug. If you really like us send us a
minimal compilable and runnable example.
* contacts:
* 	- <a href="mailto:patrick.zulian@gmail.com"> patrick.zulian@gmail.com </a> for general bugs;
* 	- <a href="mailto:alena.kopanicakova@usi.ch"> alena.kopanicakova@usi.ch </a> for solvers related bugs.
*
* We will set-up a public repository soon complete with bugtracker.
*/

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Solvers.hpp"
#include "utopia_Utils.hpp"

#ifdef WITH_TRILINOS
// Macros defined in a petsc header are conflicting with Enums defined in Trilinos
// Trilinos inclusions should always preceed petsc inclusions
#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"
#endif  // WITH_TRILINOS

#ifdef WITH_BLAS
#include "utopia_blas.hpp"
#endif  // WITH_BLAS

#ifdef WITH_PETSC
#include "utopia_petsc.hpp"
#endif  // WITH_PETSC

#ifdef WITH_CUDA
#include "utopia_cuda.hpp"
#endif  // WITH_CUDA

#ifdef WITH_UTOPIA_OPENCL
#include "utopia_opencl.hpp"
#endif  // WITH_UTOPIA_OPENCL

#ifdef WITH_M3ELINSOL
#include "utopia_m3elinsol.hpp"
#endif  // WITH_M3ELINSOL

#endif  // UTOPIA_HPP
