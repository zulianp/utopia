#ifndef UTOPIA_PETSC_DECOMPOSE_HPP
#define UTOPIA_PETSC_DECOMPOSE_HPP

#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Traits.hpp"

namespace utopia {
    bool decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions);
    bool parallel_decompose(const PetscMatrix &matrix, const int num_partitions, int *partitions);

    bool partitions_to_permutations(const PetscMatrix &matrix,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index);

    bool partitions_to_permutations(const Communicator &comm,
                                    const ArrayView<const PetscInt> &rrs,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index,
                                    std::vector<int> &rpartitions);

    bool redistribute_from_permutation(const PetscMatrix &in,
                                       const Traits<PetscMatrix>::IndexArray &permutation,
                                       PetscMatrix &out,
                                       MatReuse reuse = MAT_INITIAL_MATRIX);

    bool redistribute_from_permutation(const PetscVector &in,
                                       const Traits<PetscMatrix>::IndexArray &permutation,
                                       PetscVector &out);

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   std::vector<int> &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation);

    bool initialize_rebalance(const PetscMatrix &in,
                              std::vector<int> &partitioning,
                              Traits<PetscMatrix>::IndexArray &permutation,
                              std::vector<int> &r_partitioning,
                              Traits<PetscMatrix>::IndexArray &r_permutation);

    bool initialize_rebalance_block(
                              const int block_size,
                              const PetscMatrix &in,
                              std::vector<int> &partitioning,
                              Traits<PetscMatrix>::IndexArray &permutation,
                              std::vector<int> &r_partitioning,
                              Traits<PetscMatrix>::IndexArray &r_permutation);

    bool rebalance(const PetscMatrix &in,
                   PetscMatrix &out,
                   std::vector<int> &partitioning,
                   Traits<PetscMatrix>::IndexArray &permutation,
                   std::vector<int> &r_partitioning,
                   Traits<PetscMatrix>::IndexArray &r_permutation);

    bool partitions_to_permutations(const PetscMatrix &matrix,
                                    const int *partitions,
                                    Traits<PetscMatrix>::IndexArray &index);

    bool inverse_partition_mapping(const int comm_size,
                                   const ArrayView<const PetscInt> &original_row_ranges,
                                   const Traits<PetscMatrix>::IndexArray &permutation,
                                   std::vector<int> &partitioning);

}  // namespace utopia

#endif  // UTOPIA_PETSC_DECOMPOSE_HPP
