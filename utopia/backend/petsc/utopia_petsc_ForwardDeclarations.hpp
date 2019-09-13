#ifndef UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
#define UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include <vector>

namespace utopia {

    class PetscMatrix;
    class PetscVector;
    // class PetscIndexSet;
    //FIXME
    using PetscIndexSet = std::vector<int>;
    class PetscTraits;

    //class PetscArray
    //FIXME
    template<typename T>
    using PetscArray = std::vector<T>;

    template<class Matrix, class Vector>
    class TaoSolver;
}

#endif //UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

