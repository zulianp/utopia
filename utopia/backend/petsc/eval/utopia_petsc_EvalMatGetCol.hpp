#ifndef UTOPIA_PETSC_MAT_GET_COL_HPP
#define UTOPIA_PETSC_MAT_GET_COL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template<class Matrix, class Vector>
    class EvalGetCol<Matrix, Vector, PETSC> {
    public:
        static void apply(const Tensor<Matrix, 2> &M, Tensor<Vector, 1> &v, typename utopia::Traits<Vector>::SizeType col_id)
        {
            M.derived().col(col_id, v.derived());
        }
    };
}

#endif //UTOPIA_PETSC_MAT_GET_COL_HPP
