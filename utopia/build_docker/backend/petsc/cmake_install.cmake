# Install script for directory: /shared/utopia/backend/petsc

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
    "/shared/utopia/backend/petsc/./utopia_PetscIS.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Base.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Communicator.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Each.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Error.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_ForwardDeclarations.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_IndexSet.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Layout.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Library.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Matrix.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Matrix_impl.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Redundant.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_RowView.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Traits.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Types.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Vector.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_Vector_impl.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_debug.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_impl.hpp"
    "/shared/utopia/backend/petsc/./utopia_petsc_quirks.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Assert.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_EvalDotVecVecs.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_EvalMatGetCol.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Blocks.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Chop.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Cond.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Cond_impl.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_DotOpDot.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Factory.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Inverse.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_KroneckerProduct.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Misc.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Monitor.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Multiply.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_NZZXRow.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Parallel.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Rename.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_Residual.hpp"
    "/shared/utopia/backend/petsc/eval/utopia_petsc_Eval_VecUniqueSortSerial.hpp"
    "/shared/utopia/backend/petsc/solvers/utopia_petsc_Solver_Traits.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_ConvergedReason.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_Factorization.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_Factorizations.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_GMRES.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_KSPSolver.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_KSPSolver_impl.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_KSPSolvers.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_KSPTypes.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_LinearSolverFactory.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_RedundantLinearSolver.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_build_ksp.hpp"
    "/shared/utopia/backend/petsc/solvers/linear/utopia_petsc_solvers.hpp"
    "/shared/utopia/backend/petsc/solvers/multilevel/utopia_petsc_Multigrid.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/utopia_petsc_Newton.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/utopia_petsc_SNES.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/utopia_petsc_SNESFunction.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/utopia_petsc_SNES_impl.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/utopia_petsc_TaoSolver.hpp"
    "/shared/utopia/backend/petsc/solvers/smoothers/utopia_petsc_ConjugateGradient.hpp"
    "/shared/utopia/backend/petsc/solvers/smoothers/utopia_petsc_GaussSeidel.hpp"
    "/shared/utopia/backend/petsc/solvers/smoothers/utopia_petsc_NonLinearSmoothers.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/constrained/quadratic_programming/utopia_petsc_RedundantQPSolver.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/constrained/quadratic_programming/utopia_petsc_SemismoothNewton.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/constrained/quadratic_programming/utopia_petsc_SemismoothNewton_impl.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/constrained/quadratic_programming/utopia_petsc_TaoQPSolver.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/constrained/quadratic_programming/utopia_petsc_TaoQPSolver_impl.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/trust_region/utopia_petsc_KSPTR.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/trust_region/utopia_petsc_TaoTRQP.hpp"
    "/shared/utopia/backend/petsc/solvers/nonlinear/trust_region/utopia_petsc_TrustRegionFactory.hpp"
    )
endif()

