#include "utopia_ContactStabilizedNewmark.hpp"

#ifndef UTOPIA_ENABLE_TRILINOS_ALGEBRA

namespace utopia {
    template class ContactStabilizedNewmark<USparseMatrix, UVector>;
}

#endif  // UTOPIA_ENABLE_TRILINOS_ALGEBRA
