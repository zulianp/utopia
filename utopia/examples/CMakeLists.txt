

# Hello world
add_executable(example_hello_world EXCLUDE_FROM_ALL example_hello_world.cpp)
utopia_link_default_targets(example_hello_world)
target_link_libraries(example_hello_world PRIVATE utopia)


# Assemble and solve a linear system
add_executable(example_solve_PDE_1D EXCLUDE_FROM_ALL example_solve_PDE_1D.cpp)
utopia_link_default_targets(example_solve_PDE_1D)
target_link_libraries(example_solve_PDE_1D PRIVATE utopia)

# Loops
add_executable(example_loops EXCLUDE_FROM_ALL example_loops.cpp)
utopia_link_default_targets(example_loops)
target_link_libraries(example_loops PRIVATE utopia)

# Non-linear
add_executable(example_nonlinear_solve EXCLUDE_FROM_ALL example_nonlinear_solve.cpp)
utopia_link_default_targets(example_nonlinear_solve)
target_link_libraries(example_nonlinear_solve PRIVATE utopia)

# Communicator
add_executable(example_communicator EXCLUDE_FROM_ALL example_communicator.cpp)
utopia_link_default_targets(example_communicator)
target_link_libraries(example_communicator PRIVATE utopia)

# Simple linear system
add_executable(example_trivial_linear_system EXCLUDE_FROM_ALL example_trivial_linear_system.cpp)
utopia_link_default_targets(example_trivial_linear_system)
target_link_libraries(example_trivial_linear_system PRIVATE utopia)

# Sum of two vectors
add_executable(example_sum_of_two_vectors EXCLUDE_FROM_ALL example_sum_of_two_vectors.cpp)
utopia_link_default_targets(example_sum_of_two_vectors)
target_link_libraries(example_sum_of_two_vectors PRIVATE utopia)

# Linear solvers comparision
add_executable(example_comparison_linear_solvers EXCLUDE_FROM_ALL example_comparison_linear_solvers.cpp)
utopia_link_default_targets(example_comparison_linear_solvers)
target_link_libraries(example_comparison_linear_solvers PRIVATE utopia)

add_executable(example_linear_regression EXCLUDE_FROM_ALL example_linear_regression.cpp)
utopia_link_default_targets(example_linear_regression)
target_link_libraries(example_linear_regression PRIVATE utopia)


# Add more examples here and add target to dependecy list below
# ...


add_custom_target(utopia_examples)
add_dependencies(utopia_examples
            example_hello_world
            example_solve_PDE_1D
            example_loops
            example_nonlinear_solve
            example_communicator
            example_sum_of_two_vectors
            example_trivial_linear_system
            example_comparison_linear_solvers
            example_linear_regression
            )
