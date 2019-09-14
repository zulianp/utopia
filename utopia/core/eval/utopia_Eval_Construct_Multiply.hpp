#ifndef UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP
#define UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP


#include "utopia_Eval_Empty.hpp"
#include <cassert>

namespace utopia {

    template<class CLeft, class Left, class Right, class Traits, int Backend>
    class Eval< Construct<CLeft, Multiply<Tensor<Left, 2>, Right> >, Traits, Backend> {
    public:
        typedef utopia::Construct<CLeft, Multiply<Tensor<Left, 2>, Right> > Expr;
        typedef typename TypeAndFill<Traits, CLeft>::Type Result;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto & result      = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().derived();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

            left.multiply(right, result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template<class CLeft, class Left, class Right, class Traits, int Backend>
    class Eval< Assign<CLeft, Multiply<Tensor<Left, 2>, Right> >, Traits, Backend> {
    public:
        typedef utopia::Assign<CLeft, Multiply<Tensor<Left, 2>, Right> > Expr;
        typedef typename TypeAndFill<Traits, CLeft>::Type Result;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto & result      = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().derived();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

            left.multiply(right, result);

            UTOPIA_TRACE_END(expr);
        }
    };

}

#endif //UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP
