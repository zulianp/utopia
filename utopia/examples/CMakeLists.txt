

# Hello world
add_executable(example_hello_world EXCLUDE_FROM_ALL example_hello_world.cpp)
utopia_link_default_targets(example_hello_world)
target_link_libraries(example_hello_world utopia)


# Assemble and solve a linear system
add_executable(example_solve_PDE_1D EXCLUDE_FROM_ALL example_solve_PDE_1D.cpp)
utopia_link_default_targets(example_solve_PDE_1D)
target_link_libraries(example_solve_PDE_1D utopia)

# Loops
add_executable(example_loops EXCLUDE_FROM_ALL example_loops.cpp)
utopia_link_default_targets(example_loops)
target_link_libraries(example_loops utopia)

# Non-linear
add_executable(example_nonlinear_solve EXCLUDE_FROM_ALL example_nonlinear_solve.cpp)
utopia_link_default_targets(example_nonlinear_solve)
target_link_libraries(example_nonlinear_solve utopia)

# Communicator
add_executable(example_communicator EXCLUDE_FROM_ALL example_communicator.cpp)
utopia_link_default_targets(example_communicator)
target_link_libraries(example_communicator utopia)

# Simple linear system
add_executable(example_trivial_linear_system EXCLUDE_FROM_ALL example_trivial_linear_system.cpp)
utopia_link_default_targets(example_trivial_linear_system)
target_link_libraries(example_trivial_linear_system utopia)

# Sum of two vectors
add_executable(example_sum_of_two_vectors EXCLUDE_FROM_ALL example_sum_of_two_vectors.cpp)
utopia_link_default_targets(example_sum_of_two_vectors)
target_link_libraries(example_sum_of_two_vectors utopia)

# Add more examples here and add target to dependecy list below
# ...


add_custom_target(utopia_examples)
add_dependencies(utopia_examples
            example_hello_world
            example_solve_PDE_1D
            example_loops
            example_nonlinear_solve
            example_communicator
<<<<<<< HEAD
            example_vector_times_matrix
            example_sum_of_two_vectors
=======
            example_trivial_linear_system
>>>>>>> ee1f4f0ae349bc670ca32f2f9f12f2ffc5089cb5
            )