list(APPEND SOLVERS_MODULES
    .
    eigensolvers
    linear
    multilevel
    multilevel/transfer
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

set(LOCAL_HEADERS "")
set(LOCAL_SOURCES "")
find_project_files(CMAKE_CURRENT_SOURCE_DIR SOLVERS_MODULES LOCAL_HEADERS LOCAL_SOURCES)

foreach(MODULE ${SOLVER_MODULES})
    include_directories(${MODULE})
endforeach(MODULE)

list(APPEND LOCAL_SOURCES utopia_solvers.cpp)
add_library(utopia_solvers STATIC ${LOCAL_SOURCES})
utopia_link_default_targets(utopia_solvers)

target_include_directories(utopia_solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR})

foreach(MODULE ${SOLVERS_MODULES})
    target_include_directories(utopia_solvers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE})
endforeach(MODULE)

target_link_libraries(utopia_solvers utopia_utils utopia_core)

install(TARGETS utopia_solvers
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(FILES ${LOCAL_HEADERS} DESTINATION include)

set(UTOPIA_COMPONENTS "${UTOPIA_COMPONENTS};utopia_solvers" PARENT_SCOPE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UTOPIA_DEV_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UTOPIA_SPECIAL_DEBUG_FLAGS} ${UTOPIA_DEV_FLAGS}")
