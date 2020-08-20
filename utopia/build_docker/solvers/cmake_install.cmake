# Install script for directory: /shared/utopia/solvers

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
    "/shared/utopia/solvers/./utopia_ConvergenceReason.hpp"
    "/shared/utopia/solvers/./utopia_Memory.hpp"
    "/shared/utopia/solvers/./utopia_Monitor.hpp"
    "/shared/utopia/solvers/./utopia_NumericalTollerance.hpp"
    "/shared/utopia/solvers/./utopia_PrintInfo.hpp"
    "/shared/utopia/solvers/./utopia_SolutionStatus.hpp"
    "/shared/utopia/solvers/./utopia_SolverForwardDeclarations.hpp"
    "/shared/utopia/solvers/./utopia_SolverType.hpp"
    "/shared/utopia/solvers/./utopia_Solvers.hpp"
    "/shared/utopia/solvers/eigensolvers/utopia_EigenSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_AbstractLinearSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_CrossBackendLinearSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_DirectSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_IterativeSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_Linear.hpp"
    "/shared/utopia/solvers/linear/utopia_LinearSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_LinearSolverFactory.hpp"
    "/shared/utopia/solvers/linear/utopia_LinearSolverInterfaces.hpp"
    "/shared/utopia/solvers/linear/utopia_MatrixFreeLinearSolver.hpp"
    "/shared/utopia/solvers/linear/utopia_PreconditionedSolver.hpp"
    "/shared/utopia/solvers/multilevel/utopia_FAS.hpp"
    "/shared/utopia/solvers/multilevel/utopia_Level.hpp"
    "/shared/utopia/solvers/multilevel/utopia_LevelMemory.hpp"
    "/shared/utopia/solvers/multilevel/utopia_LinearMultiLevel.hpp"
    "/shared/utopia/solvers/multilevel/utopia_MG_OPT.hpp"
    "/shared/utopia/solvers/multilevel/utopia_MultiLevelBase.hpp"
    "/shared/utopia/solvers/multilevel/utopia_MultiLevelMask.hpp"
    "/shared/utopia/solvers/multilevel/utopia_Multigrid.hpp"
    "/shared/utopia/solvers/multilevel/utopia_MultigridQR.hpp"
    "/shared/utopia/solvers/multilevel/utopia_Multilevel.hpp"
    "/shared/utopia/solvers/multilevel/utopia_NonlinearMultiLevelBase.hpp"
    "/shared/utopia/solvers/multilevel/utopia_QuasiRMTR.hpp"
    "/shared/utopia/solvers/multilevel/utopia_QuasiRMTR_inf.hpp"
    "/shared/utopia/solvers/multilevel/utopia_RMTR.hpp"
    "/shared/utopia/solvers/multilevel/utopia_RMTRBase.hpp"
    "/shared/utopia/solvers/multilevel/utopia_RMTRParams.hpp"
    "/shared/utopia/solvers/multilevel/utopia_RMTRVcycleImpl.hpp"
    "/shared/utopia/solvers/multilevel/utopia_RMTR_inf.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_IPTransfer.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_IPTransferNested.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_IdentityTransfer.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_MatrixTransfer.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_MatrixTruncatedTransfer.hpp"
    "/shared/utopia/solvers/multilevel/transfer/utopia_Transfer.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_BoxGelmanMandel.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_BoxKornhuber.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_IdentityConstraints.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_MLConstraintsIncludes.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_MultiLevelVariableBoundInterface.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_TRBoundsGelmanMandel.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_TRBoundsGratton.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_TRBoundsKornhuber.hpp"
    "/shared/utopia/solvers/multilevel/constraints/utopia_TRBoxMixConstraints.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_FunEvalsIncludes.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MLEvalFirstOrder.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MLEvalFirstOrderDF.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MLEvalFirstOrderMGOPT.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MLEvalGalerkin.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MLEvalSecondOrder.hpp"
    "/shared/utopia/solvers/multilevel/fun_evals/utopia_MultiLevelEvaluations.hpp"
    "/shared/utopia/solvers/smoothers/utopia_BiCGStab.hpp"
    "/shared/utopia/solvers/smoothers/utopia_BiCGStab_impl.hpp"
    "/shared/utopia/solvers/smoothers/utopia_ConjugateGradient.hpp"
    "/shared/utopia/solvers/smoothers/utopia_ConjugateGradient_impl.hpp"
    "/shared/utopia/solvers/smoothers/utopia_GaussSeidel.hpp"
    "/shared/utopia/solvers/smoothers/utopia_Jacobi.hpp"
    "/shared/utopia/solvers/smoothers/utopia_NonLinearJacobi.hpp"
    "/shared/utopia/solvers/smoothers/utopia_NonLinearSmoother.hpp"
    "/shared/utopia/solvers/smoothers/utopia_PointJacobi.hpp"
    "/shared/utopia/solvers/smoothers/utopia_Smoother.hpp"
    "/shared/utopia/solvers/smoothers/utopia_Smoothers.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_ExtendedFunction.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_Function.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_FunctionNormalEq.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_GradientDescent.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_MSSolver.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_Newton.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_NewtonBase.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_NonLinearSolver.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_Nonlinear.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_NonlinearLeastSquaresSolver.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_NonlinearSolverFactory.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_NonlinearSolverInterfaces.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_QuadraticFunction.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_QuasiNewton.hpp"
    "/shared/utopia/solvers/nonlinear/utopia_QuasiNewtonBase.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_ASTRUM.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_AffineSimilarity.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_AffineSimilarityAW.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_LevenbergMarquardt.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_PseudoContinuation.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_PseudoContinuationIncludes.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_PseudoTrustRegion.hpp"
    "/shared/utopia/solvers/nonlinear/pseudo_transient_continuation/utopia_RosenbrockTrustRegion.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/utopia_BoxConstraints.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/utopia_ConstrainedIncludes.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/utopia_NonlinSemismoothNewton.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/utopia_QuasiNewtonBound.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/utopia_VariableBoundSolverInterface.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_BlockQPSolver.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_BlockQPSolver_impl.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_GenericSemismoothNewton.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_MPRGP.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedConjugateGradient.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGaussSeidel.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGaussSeidelNew.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGaussSeidelQR.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGaussSeidelSweep.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGaussSeidel_impl.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_ProjectedGradient.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_QPSolver.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_RedundantQPSolver.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_SemismoothNewton.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_SemismoothNewton_impl.hpp"
    "/shared/utopia/solvers/nonlinear/constrained/quadratic_programming/utopia_SemismoothNewton_old.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_Backtracking.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_LS_Factory.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_LS_Strategy.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_LS_normal_eq.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_LineSearchIncludes.hpp"
    "/shared/utopia/solvers/nonlinear/line_search/utopia_SimpleBacktracking.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_QuasiTrustRegion.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_QuasiTrustRegionVariableBound.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TRBase.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TRNormalEquation.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TrustRegion.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TrustRegionFactory.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TrustRegionIncludes.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/utopia_TrustRegionVariableBound.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_CauchyPoint.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_Dogleg.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_MoreSorensen.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_SteihaugToint.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_TRQuadraticFunction.hpp"
    "/shared/utopia/solvers/nonlinear/trust_region/QP_subproblems/utopia_TRSubproblem.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_BFGS.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_HessianApproximation.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_HessianApproximations.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_JFNK.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_LBFGS.hpp"
    "/shared/utopia/solvers/nonlinear/hessian_approximations/utopia_LSR1.hpp"
    "/shared/utopia/solvers/preconditioners/utopia_BlockPreconditioner.hpp"
    "/shared/utopia/solvers/preconditioners/utopia_Preconditioner.hpp"
    "/shared/utopia/solvers/preconditioners/utopia_Preconditioners.hpp"
    "/shared/utopia/solvers/preconditioners/utopia_TrivialPreconditioners.hpp"
    "/shared/utopia/solvers/saddlepoint/utopia_SPBlockConjugateGradient.hpp"
    "/shared/utopia/solvers/saddlepoint/utopia_SPStaticCondensation.hpp"
    "/shared/utopia/solvers/saddlepoint/utopia_SPStaticCondensationKrylov.hpp"
    "/shared/utopia/solvers/saddlepoint/utopia_SaddlePoint.hpp"
    )
endif()

