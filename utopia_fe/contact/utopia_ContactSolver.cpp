#include "utopia_ContactSolver.hpp"

//FIXME
#ifndef WITH_TRILINOS_ALGEBRA
namespace utopia {
    template class ContactSolver<USparseMatrix, UVector>;
}
#endif //WITH_TRILINOS_ALGEBRA
