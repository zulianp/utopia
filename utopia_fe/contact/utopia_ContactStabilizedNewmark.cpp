#include "utopia_ContactStabilizedNewmark.hpp"

#ifndef WITH_TRILINOS_ALGEBRA

namespace utopia {
    template class ContactStabilizedNewmark<USparseMatrix, UVector>;
}

#endif  // WITH_TRILINOS_ALGEBRA
