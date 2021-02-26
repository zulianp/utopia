#ifndef UTOPIA_PETSC_DILU_ALGORITHM_HPP
#define UTOPIA_PETSC_DILU_ALGORITHM_HPP

#include <memory>

#include "utopia_petsc_Types.hpp"

#include "utopia_DILUDecompose.hpp"

namespace utopia {

    template <>
    class DILUAlgorithm<PetscMatrix, PetscVector> final : public ILUAlgorithm<PetscMatrix, PetscVector> {
    public:
        bool update(const PetscMatrix &mat) override;
        void apply(const PetscVector &b, PetscVector &x) override;
        void read(Input &) override;

        DILUAlgorithm();
        ~DILUAlgorithm();

        class Impl;

    private:
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_DILU_ALGORITHM_HPP
