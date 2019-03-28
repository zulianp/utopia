#ifndef UTOPIA_PETSC_MAT_GET_COL_HPP
#define UTOPIA_PETSC_MAT_GET_COL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia
{

    template<class Matrix, class Vector>
    class EvalGetCol<Matrix, Vector, PETSC>
    {
        public:
            static void apply(const Wrapper<Matrix, 2> &M, Wrapper<Vector, 1> &v, typename utopia::Traits<Vector>::SizeType col_id)
            {
                Backend<typename Traits<Matrix>::Scalar, Traits<Matrix>::Backend>::Instance().mat_get_col(M.implementation(), v.implementation(), col_id);
            }
    };


}


#endif //UTOPIA_PETSC_MAT_GET_COL_HPP
