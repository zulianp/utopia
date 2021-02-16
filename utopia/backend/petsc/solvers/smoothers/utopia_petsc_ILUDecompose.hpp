#ifndef UTOPIA_PETSC_ILU_DECOMPOSE_HPP
#define UTOPIA_PETSC_ILU_DECOMPOSE_HPP

#include "utopia_petsc_Types.hpp"

#include "utopia_CRSMatrix.hpp"
#include "utopia_ILU.hpp"

namespace utopia {

    template <int BlockSize>
    void crs_block_matrix(const PetscMatrix &in,
                          CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out);

    template <int BlockSize>
    void crs_block_matrix_split_diag(const PetscMatrix &in,
                                     CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out,
                                     std::vector<PetscReal> &diag);

    template <int BlockSize>
    void crs_block_matrix_update(const PetscMatrix &in,
                                 CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> &out,
                                 std::vector<PetscReal> &diag);

    template <>
    class ILUDecompose<PetscMatrix, PETSC> final : public ILUAlgorithm<PetscMatrix, PetscVector> {
    public:
        static bool decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified);
        static void apply(const PetscMatrix &ilu, const PetscVector &b, PetscVector &x);

        bool update(const PetscMatrix &mat) override;
        void apply(const PetscVector &b, PetscVector &x) override;
        void read(Input &) override;

    private:
        PetscMatrix ilu_;
        bool modified_{false};
    };

    template <int BlockSize>
    class BlockILUAlgorithm<PetscMatrix, BlockSize> final : public ILUAlgorithm<PetscMatrix, PetscVector> {
    public:
        bool update(const PetscMatrix &mat) override;
        void apply(const PetscVector &b, PetscVector &x) override;

        void read(Input &) override;

    private:
        CRSMatrix<std::vector<PetscScalar>, std::vector<PetscInt>, BlockSize> ilu_;
        PetscVector L_inv_b_;
        bool print_matrices_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_ILU_DECOMPOSE_HPP
