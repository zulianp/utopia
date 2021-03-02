#ifndef UTOPIA_PETSC_DILU_ALGORITHM_HPP
#define UTOPIA_PETSC_DILU_ALGORITHM_HPP

#include <memory>

#include "utopia_petsc_Types.hpp"

#include "utopia_DILUDecompose.hpp"

namespace utopia {

    template <>
    class DILUAlgorithm<PetscMatrix, PetscVector> final : public ILUAlgorithm<PetscMatrix, PetscVector> {
    public:
        using Super = utopia::ILUAlgorithm<PetscMatrix, PetscVector>;

        bool update(const PetscMatrix &mat) override;
        void apply(const PetscVector &b, PetscVector &x) override;
        void read(Input &) override;

        inline DILUAlgorithm *clone() const override { return new DILUAlgorithm(); }

        DILUAlgorithm();
        ~DILUAlgorithm();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        void init();
    };

    template <int BlockSize>
    class BlockDILUAlgorithm<PetscMatrix, PetscVector, BlockSize> final
        : public ILUAlgorithm<PetscMatrix, PetscVector> {
    public:
        using Super = utopia::ILUAlgorithm<PetscMatrix, PetscVector>;

        bool update(const PetscMatrix &mat) override;
        void apply(const PetscVector &b, PetscVector &x) override;
        void read(Input &) override;

        inline BlockDILUAlgorithm *clone() const override { return new BlockDILUAlgorithm(); }

        BlockDILUAlgorithm();
        ~BlockDILUAlgorithm();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        void init();
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_DILU_ALGORITHM_HPP
