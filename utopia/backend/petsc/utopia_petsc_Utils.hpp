#ifndef UTOPIA_PETSC_UTILS_HPP
#define UTOPIA_PETSC_UTILS_HPP

#include "utopia_petsc_Types.hpp"

namespace utopia {
    void optimize_nnz(PetscMatrix &A);
    bool is_diagonally_dominant(const PetscMatrix &A);
    void local_block_view(const PetscMatrix &mat, PetscMatrix &block);
    // void local_view(const PetscVector &vec, PetscVector &lv);
}  // namespace utopia

#endif  // UTOPIA_PETSC_UTILS_HPP
