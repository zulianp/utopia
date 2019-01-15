#ifndef UTOPIA_PETSC_EVAL_BLOCKS_HPP
#define UTOPIA_PETSC_EVAL_BLOCKS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Blocks.hpp"

namespace utopia {

    template<class Left, class Right, class Traits>
    class Eval< Construct<Wrapper<Left, 1>, Blocks<Right> >, Traits, PETSC> {
    public:
        inline static bool apply(const Construct<Wrapper<Left, 1>, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).build_blocks(
                Eval<Wrapper<Left, 1>,  Traits>::apply(expr.left()),
                expr.right()
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits>
    class Eval< Construct<Wrapper<Left, 2>, Blocks<Right> >, Traits, PETSC> {
    public:
        inline static bool apply(const Construct<Wrapper<Left, 2>, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            if(is_sparse<Right>::value) {
                UTOPIA_BACKEND(Traits).build_blocks(
                    Eval<Wrapper<Left, 2>,  Traits>::apply(expr.left()),
                    expr.right()
                );
            } else {
                Eval<Construct<Wrapper<Left, 2>, Blocks<Right> >, Traits, HOMEMADE>::apply(expr);
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_PETSC_EVAL_BLOCKS_HPP
