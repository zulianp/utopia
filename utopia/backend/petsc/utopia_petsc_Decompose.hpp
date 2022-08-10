#ifndef UTOPIA_PETSC_DECOMPOSE_HPP
#define UTOPIA_PETSC_DECOMPOSE_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {
    bool decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions);
    bool parallel_decompose(const PetscMatrix &matrix, const int numparitions, int *paritions);
}

#endif  // UTOPIA_PETSC_DECOMPOSE_HPP
