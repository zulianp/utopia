#ifndef UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
#define UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_petsc_Base.hpp"

#include <vector>

#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 11, 3)
#include <petscsystypes.h>
#else
// FIXME find the correct header
#include "petscvec.h"
#endif

namespace utopia {

    class PetscMatrix;
    class PetscVector;
    // class PetscIndexSet;
    // FIXME
    using PetscIndexSet = std::vector<PetscInt>;
    class PetscTraits;

    // class PetscArray
    // FIXME
    template <typename T>
    using PetscArray = std::vector<T>;

    template <class Matrix, class Vector>
    class TaoSolver;

    template <typename Matrix, typename Vector, int Backend>
    class KSPSolver;

    class PetscCommunicator;
}  // namespace utopia

#endif  // UTOPIA_PETSC_FORWARD_DECLARATIONS_HPP
