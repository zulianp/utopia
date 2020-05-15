//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_NUMBER_HPP
#define UTOPIA_UTOPIA_EVAL_NUMBER_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Literal.hpp"

namespace utopia {

    template <typename L, class Traits, int Backend>
    class Eval<Number<L>, Traits, Backend> {
    public:
        inline static L apply(const Number<L> &num) { return num; }
    };

    template <class Expr, class Traits, int Backend>
    class Eval<Boolean<Expr>, Traits, Backend> {
    public:
        inline static bool apply(const Boolean<Expr> &expr) { return Eval<Expr, Traits>(expr.expr()); }
    };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_NUMBER_HPP
