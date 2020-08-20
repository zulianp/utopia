# Install script for directory: /shared/utopia/backend/blas

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/shared/utopia/backend/blas/./utopia_blas.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Algorithms.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Array.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_ForwardDeclarations.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_IndexSet.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Matrix.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_RowView.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Traits.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_TrustRegionFactory.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Types.hpp"
    "/shared/utopia/backend/blas/./utopia_blas_Vector.hpp"
    "/shared/utopia/backend/blas/solvers/utopia_blas_LinearSolverFactory.hpp"
    "/shared/utopia/backend/blas/solvers/utopia_blas_solvers.hpp"
    "/shared/utopia/backend/blas/eigensolvers/utopia_blas_eigensolvers.hpp"
    "/shared/utopia/backend/blas/eval/utopia_blas_Each.hpp"
    "/shared/utopia/backend/blas/eval/utopia_blas_Eval.hpp"
    "/shared/utopia/backend/blas/eval/utopia_blas_Eval_KroneckerProduct.hpp"
    "/shared/utopia/backend/blas/eval/utopia_blas_Eval_Misc.hpp"
    "/shared/utopia/backend/blas/eval/utopia_blas_Eval_Multiply.hpp"
    "/shared/utopia/backend/blas/solvers/linear/lapack/utopia_Lapack.hpp"
    "/shared/utopia/backend/blas/eigensolvers/lapack/utopia_LapackEigenSolver.hpp"
    "/shared/utopia/backend/blas/eigensolvers/lapack/utopia_lapack_Algorithms.hpp"
    )
endif()

