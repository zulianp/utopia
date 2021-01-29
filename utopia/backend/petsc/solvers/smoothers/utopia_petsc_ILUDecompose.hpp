#ifndef UTOPIA_PETSC_ILU_DECOMPOSE_HPP
#define UTOPIA_PETSC_ILU_DECOMPOSE_HPP

#include "utopia_petsc_Types.hpp"

#include "utopia_ILU.hpp"

namespace utopia {

    template <>
    class ILUDecompose<PetscMatrix, PETSC> {
    public:
        static void decompose(const PetscMatrix &mat, PetscMatrix &out, const bool modified);
        static void apply(const PetscMatrix &ilu, const PetscVector &b, PetscVector &x);

        // static void apply_vi(const PetscMatrix &ilu,
        //                      const PetscVector &lb,
        //                      const PetscVector &ub,
        //                      const PetscVector &b,
        //                      PetscVector &x);
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_ILU_DECOMPOSE_HPP
