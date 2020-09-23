#ifndef UTOPIA_PETSC_COND_HPP
#define UTOPIA_PETSC_COND_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_SLEPC

#include "utopia_Cond.hpp"

namespace utopia {

    template <class Matrix>
    class Cond<Matrix, PETSC> {
    public:
        using Scalar = UTOPIA_SCALAR(Matrix);
        static Scalar apply(const Matrix &H);
    };
}  // namespace utopia

#endif  // UTOPIA_WITH_SLEPC
#endif  // UTOPIA_PETSC_COND_HPP
