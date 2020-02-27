#ifndef UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
#define UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include <vector>
#include <petscsystypes.h>

namespace utopia {

    class PetscMatrix;
    class PetscVector;
    // class PetscIndexSet;
    //FIXME
    using PetscIndexSet = std::vector<PetscInt>;
    class PetscTraits;

    //class PetscArray
    //FIXME
    template<typename T>
    using PetscArray = std::vector<T>;

    template<class Matrix, class Vector>
    class TaoSolver;

    template<typename Matrix, typename Vector, int Backend>
    class KSPSolver;

    class PetscCommunicator;
}

#endif //UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

