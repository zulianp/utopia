#include "utopia_ContactStabilizedNewmark.hpp"

#ifndef UTOPIA_WITH_TRILINOS_ALGEBRA

namespace utopia {
    template class ContactStabilizedNewmark<USparseMatrix, UVector>;
}

#endif  // UTOPIA_WITH_TRILINOS_ALGEBRA
