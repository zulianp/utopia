#ifndef UTOPIA_LINEAR_SOLVE_EXPORT_HPP
#define UTOPIA_LINEAR_SOLVE_EXPORT_HPP

typedef void * USolver;
typedef void * UMat;
typedef void * UVec;

void utopia_initialize(int argc, char *argv[]);
void utopia_finalize();
void utopia_create_solver(const char name[], USolver * solver);
void utopia_destroy_solver(USolver * solver);
void utopia_solve(USolver solver, UMat A, UVec b, UVec x);
void utopia_print_solver_info(USolver ptr);
void utopia_print_version();

#endif
