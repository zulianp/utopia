
#ifndef UTOPIA_UTOPIA_PETSCERROR_HPP
#define UTOPIA_UTOPIA_PETSCERROR_HPP

#include <iostream>
#include "utopia_Base.hpp"
#include "utopia_IOStream.hpp"

#include <petscsys.h>

namespace utopia {
    class PetscErrorHandler {
    public:
        static bool Check(PetscErrorCode err) {
            if (err) utopia::out() << "PetscErrorCode=" << err << std::endl;
            return true;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_PETSCERROR_HPP
