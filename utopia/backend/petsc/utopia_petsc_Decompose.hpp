#ifndef UTOPIA_PETSC_DECOMPOSE_HPP
#define UTOPIA_PETSC_DECOMPOSE_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Traits.hpp"

namespace utopia {
    bool decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions);
    bool parallel_decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions);

    // bool partitions_to_permutations(const PetscMatrix &matrix,
    //                                 const int *partitions,
    //                                 Traits<PetscMatrix>::SizeType *index);

    bool partitions_to_permutations(const PetscMatrix &matrix,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index);

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   Traits<PetscMatrix>::IndexArray &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation);

}  // namespace utopia

#endif  // UTOPIA_PETSC_DECOMPOSE_HPP
