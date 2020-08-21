#ifndef UTOPIA_CONFIG_HPP
#define UTOPIA_CONFIG_HPP

#define UTOPIA_WITH_BLAS
#define UTOPIA_WITH_LAPACK
#define UTOPIA_WITH_PETSC
/* #undef UTOPIA_WITH_SLEPC */
/* #undef UTOPIA_WITH_TRILINOS */
/* #undef UTOPIA_WITH_TRILINOS_AMESOS2 */
/* #undef UTOPIA_WITH_TRILINOS_BELOS */
/* #undef UTOPIA_WITH_TRILINOS_IFPACK2 */
/* #undef UTOPIA_WITH_TRILINOS_MUELU */
/* #undef UTOPIA_WITH_TRILINOS_TPETRAEXT */
#define UTOPIA_WITH_CPP11
/* #undef UTOPIA_WITH_UMFPACK */
#define UTOPIA_WITH_MPI
/* #undef UTOPIA_WITH_OPEN_BLAS */
/* #undef UTOPIA_WITH_PETSC_CRAY */
/* #undef UTOPIA_WITH_CUDA */
/* #undef WITH_UTOPIA_OPENCL */
/* #undef UTOPIA_WITH_EIGEN_3 */
/* #undef UTOPIA_TRACE_ENABLED */
/* #undef UTOPIA_TRACE_EXPR_ENABLED */
/* #undef UTOPIA_ENABLE_NO_ALLOC_REGIONS*/
/* #undef WITH_PASSO_EXTENSIONS */
/* #undef UTOPIA_WITH_TINY_EXPR */
/* #undef UTOPIA_WITH_JSON */
/* #undef UTOPIA_WITH_M3ELINSOL */
#define UTOPIA_WITH_CPP14
/* #undef UTOPIA_WITH_CODE_COVERAGE */
#define UTOPIA_DEPRECATED_API

/* #undef UTOPIA_PETSC_VERSION_MAJOR */
/* #undef UTOPIA_PETSC_VERSION_MINOR */
/* #undef UTOPIA_PETSC_VERSION_SUBMINOR */

#ifdef UTOPIA_WITH_TRILINOS
#include <Kokkos_Macros.hpp>
#define UTOPIA_LAMBDA KOKKOS_LAMBDA
#define UTOPIA_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#define UTOPIA_FUNCTION KOKKOS_FUNCTION
#define NVCC_PRIVATE public:
#else
#define UTOPIA_LAMBDA [=]
#define UTOPIA_INLINE_FUNCTION inline
#define UTOPIA_FUNCTION
#define NVCC_PRIVATE
#endif

#ifdef UTOPIA_WITH_OPEN_BLAS
#define UTOPIA_WITH_BLAS
#endif

#ifdef UTOPIA_WITH_PETSC_CRAY
#define PETSC_CXX_STATIC_INLINE static inline
#define PETSC_CXX_RESTRICT __restrict__
#endif

#ifdef UTOPIA_WITH_CPP11  // Global enabler
/* #undef HAS_CXX11_AUTO */
/* #undef HAS_CXX11_NULLPTR */
/* #undef HAS_CXX11_LAMBDA */
/* #undef HAS_CXX11_STATIC_ASSERT */
/* #undef HAS_CXX11_RVALUE_REFERENCES */
/* #undef HAS_CXX11_DECLTYPE */
/* #undef HAS_CXX11_CSTDINT_H */
/* #undef HAS_CXX11_LONG_LONG */
/* #undef HAS_CXX11_VARIADIC_TEMPLATES */
/* #undef HAS_CXX11_CONSTEXPR */
/* #undef HAS_CXX11_SIZEOF_MEMBER */
/* #undef HAS_CXX11_FUNC */
#endif  // UTOPIA_WITH_CPP11

#ifndef NDEBUG
#include <assert.h>
/* #undef ENABLE_LOCK_CHECK */
#ifdef ENABLE_LOCK_CHECK
#define assert_enabled(expr_) assert(expr_)
#else
#define assert_enabled(expr_)
#endif  // ENABLE_LOCK_CHECK
#else
#define assert_enabled(expr_)
#endif  // NDEBUG

#include "utopia_compiler_features.hpp"

#endif  // UTOPIA_CONFIG_HPP
