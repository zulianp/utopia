#ifndef UTOPIA_EVAL_CHOP_HPP
#define UTOPIA_EVAL_CHOP_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Backend.hpp"
#include "utopia_Temp.hpp"

#include <cassert>

/*! @file
* Petsc language extensions
*/

namespace utopia {
    template<class Matrix>
    class ChopSmallerThan<Matrix, PETSC> {
    public:
        using Scalar   = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        static void apply(const Tensor<Matrix, 2> &A, const Scalar &eps);
    };

    template<class Matrix>
    class ChopGreaterThan<Matrix, PETSC> {
    public:
        using Scalar   = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        static void apply(const Tensor<Matrix, 2> &A, const Scalar &eps);
    };

}

#endif //UTOPIA_EVAL_CHOP_HPP
