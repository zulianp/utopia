list(APPEND SOLVERS_MODULES
    .
    eigensolvers
    linear
    multilevel
    multilevel/transfer
    multilevel/constraints
    multilevel/fun_evals
    smoothers
    nonlinear
    nonlinear/pseudo_transient_continuation
    nonlinear/constrained
    nonlinear/constrained/quadratic_programming
    nonlinear/line_search
    nonlinear/trust_region
    nonlinear/trust_region/QP_subproblems
    nonlinear/hessian_approximations
    preconditioners
    saddlepoint
)
utopia_add_library(utopia_solvers ${CMAKE_CURRENT_SOURCE_DIR} "${SOLVERS_MODULES}" "utopia_utils;utopia_core")
