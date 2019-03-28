#ifndef UTOPIA_UTOPIA_EVAL_STRUCTURE_HPP
#define UTOPIA_UTOPIA_EVAL_STRUCTURE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Structure.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Left, Structure<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Structure<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).build_from_structure(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_STRUCTURE_HPP
