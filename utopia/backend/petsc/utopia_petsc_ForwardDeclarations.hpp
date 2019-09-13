#ifndef UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
#define UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

#include <vector>

namespace utopia {

    class PetscMatrix;
    class PetscVector;
    // class PetscIndexSet;
    //FIXME
    using PetscIndexSet = std::vector<int>;
    class PetscTraits;

    template<typename T>
    class PetscArray;

    template<class Matrix, class Vector>
    class TaoSolver;
}

#endif //UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

