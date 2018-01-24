
#ifndef UTOPIA_UTOPIA_PETSCERROR_HPP
#define UTOPIA_UTOPIA_PETSCERROR_HPP

#include "utopia_Base.hpp"
#include <iostream>

#include <petscsys.h>

namespace utopia {
    class PetscErrorHandler {
    public:
        static bool Check(PetscErrorCode err) {
            if (err) std::cout << "PetscErrorCode=" << err << std::endl;
            return true;
        }
    };
}

#endif //UTOPIA_UTOPIA_PETSCERROR_HPP
