#ifndef UTOPIA_UTOPIA_EVAL_BLOCKS_HPP
#define UTOPIA_UTOPIA_EVAL_BLOCKS_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Blocks.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Left, Blocks<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Blocks<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).build_blocks(
                    Eval<Left,  Traits>::apply(expr.left()),
                    expr.right()
            );

			UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_BLOCKS_HPP
