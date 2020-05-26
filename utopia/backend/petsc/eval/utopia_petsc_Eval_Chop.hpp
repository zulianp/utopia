#ifndef UTOPIA_PETSC_MAT_CHOP_HPP
#define UTOPIA_PETSC_MAT_CHOP_HPP

#include "utopia_MatChop.hpp"
#include "utopia_Traits.hpp"
#include "utopia_petsc_Matrix.hpp"

namespace utopia {

    template <>
    class Chop<PetscMatrix, PETSC> {
    public:
        using Scalar = typename utopia::Traits<PetscMatrix>::Scalar;
        using SizeType = typename utopia::Traits<PetscMatrix>::SizeType;

        static void apply(PetscMatrix &mat, const Scalar &eps);
    };

    template <>
    class ChopSmallerThan<PetscMatrix, PETSC> {
    public:
        using Scalar = typename utopia::Traits<PetscMatrix>::Scalar;
        using SizeType = typename utopia::Traits<PetscMatrix>::SizeType;

        static void apply(PetscMatrix &mat, const Scalar &eps);
    };

    template <>
    class ChopGreaterThan<PetscMatrix, PETSC> {
    public:
        using Scalar = typename utopia::Traits<PetscMatrix>::Scalar;
        using SizeType = typename utopia::Traits<PetscMatrix>::SizeType;

        static void apply(PetscMatrix &mat, const Scalar &eps);
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_MAT_CHOP_HPP
