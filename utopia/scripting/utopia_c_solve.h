#ifndef UTOPIA_LINEAR_SOLVE_EXPORT_HPP
#define UTOPIA_LINEAR_SOLVE_EXPORT_HPP

void utopia_initialize(int argc, char *argv[]);
void utopia_finalize();
void utopia_create_solver(const char name[], void **solver);
void utopia_destroy_solver(void **solver);
void utopia_solve(void *solver, void*A, void*b, void *x);
void utopia_print_solver_info(void *ptr);
void utopia_print_version();

#endif
