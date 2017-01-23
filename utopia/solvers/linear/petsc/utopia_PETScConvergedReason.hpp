#ifndef UTOPIA_PETSC_CONVERGED_REASON_HPP
#define UTOPIA_PETSC_CONVERGED_REASON_HPP 

#include <ostream>
#include <iostream>

namespace utopia {
	void print_ksp_converged_reason(const int code, std::ostream &os = std::cout);
	bool ksp_convergence_fatal_failure(const int code);
}

#endif //UTOPIA_PETSC_CONVERGED_REASON_HPP
