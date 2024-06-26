#include "utopia_ContactSolver.hpp"

// FIXME
#ifndef UTOPIA_ENABLE_TRILINOS_ALGEBRA
namespace utopia {
    template class ContactSolver<USparseMatrix, UVector>;
}
#endif  // UTOPIA_ENABLE_TRILINOS_ALGEBRA
