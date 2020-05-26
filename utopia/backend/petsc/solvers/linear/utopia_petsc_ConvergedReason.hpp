#ifndef UTOPIA_PETSC_CONVERGED_REASON_HPP
#define UTOPIA_PETSC_CONVERGED_REASON_HPP

#include <iostream>
#include <ostream>

namespace utopia {
    void print_ksp_converged_reason(const int code, std::ostream &os = std::cout);
    bool ksp_convergence_fatal_failure(const int code);
}  // namespace utopia

#endif  // UTOPIA_PETSC_CONVERGED_REASON_HPP
