list(
    APPEND
    SOLVERS_MODULES
    .
    eigensolvers
    linear
    multilevel
    multilevel/transfer
    multilevel/constraints
    multilevel/fun_evals
    multilevel/agglomerate
    smoothers
    nonlinear
    nonlinear/pseudo_transient_continuation
    nonlinear/constrained
    nonlinear/constrained/barrier
    nonlinear/constrained/interior_point
    nonlinear/constrained/smoothers
    nonlinear/constrained/quadratic_programming
    nonlinear/line_search
    nonlinear/trust_region
    nonlinear/trust_region/QP_subproblems
    nonlinear/hessian_approximations
    preconditioners
    saddlepoint)

utopia_add_library(${CMAKE_CURRENT_SOURCE_DIR} "${SOLVERS_MODULES}")
