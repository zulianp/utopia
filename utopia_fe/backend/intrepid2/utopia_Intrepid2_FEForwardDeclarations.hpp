#ifndef UTOPIA_INTREPID2_FE_FORWARD_DECLARATIONS_HPP
#define UTOPIA_INTREPID2_FE_FORWARD_DECLARATIONS_HPP

#include "utopia_FEForwardDeclarations.hpp"

#ifdef NDEBUG
#if defined(DEBUG)
#error "DEBUG macro defined when it should not be defined"
#endif
#endif  // NDEBUG

static const int INTREPID2_TAG = 1002;

namespace utopia {
    class Intrepid2FunctionSpace;
}

#endif  // UTOPIA_INTREPID2_FE_FORWARD_DECLARATIONS_HPP