#ifndef UTOPIA_PETSC_EVAL_BLOCKS_HPP
#define UTOPIA_PETSC_EVAL_BLOCKS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Blocks.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    void build_blocks(PetscMatrix &left, const Blocks<PetscMatrix> &blocks);
    void build_blocks(PetscVector &left, const Blocks<PetscVector> &blocks);

    template<class Left, class Right, class Traits>
    class Eval< Construct<Tensor<Left, 1>, Blocks<Right> >, Traits, PETSC> {
    public:
        inline static bool apply(const Construct<Tensor<Left, 1>, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            build_blocks(Eval<Tensor<Left, 1>,  Traits>::apply(expr.left()), expr.right());

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits>
    class Eval< Construct<Tensor<Left, 2>, Blocks<Right> >, Traits, PETSC> {
    public:
        inline static bool apply(const Construct<Tensor<Left, 2>, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            build_blocks(Eval<Tensor<Left, 2>,  Traits>::apply(expr.left()),   expr.right());

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_PETSC_EVAL_BLOCKS_HPP
