#ifndef UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
#define UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP

#include "utopia_Eval_Binary.hpp"
#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Right>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Left, Right>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            EvalBinaryAux<Result>::apply(
                Eval<Left, Traits>::apply(expr.left()), Eval<Right, Traits>::apply(expr.right()), Multiplies(), result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Transposed<Left>, Right>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Transposed<Left>, Right>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&l = Eval<Left, Traits>::apply(expr.left().expr());
            auto &&r = Eval<Right, Traits>::apply(expr.right());

            l.transpose_multiply(r, result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Transposed<Right>>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Left, Transposed<Right>>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left())
                .multiply_transpose(Eval<Right, Traits>::apply(expr.right().expr()), result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Transposed<Left>, Transposed<Right>>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Transposed<Left>, Transposed<Right>>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left().expr())
                .multiply(true, true, Eval<Right, Traits>::apply(expr.right().expr()), result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Transposed<Multiply<Left, Right>>, Traits, Backend> {
    public:
        using Expr = utopia::Transposed<Multiply<Left, Right>>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.expr().right())
                .multiply(true, true, Eval<Left, Traits>::apply(expr.expr().left()), result);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Result, class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Result, Multiply<Left, Right>>, Traits, Backend> {
    public:
        using EvalMultiply = utopia::Eval<Multiply<Left, Right>, Traits, Backend>;
        using Expr = utopia::Assign<Result, Multiply<Left, Right>>;

        inline static void apply(const Expr &expr) { EvalMultiply::apply(expr.right(), expr.left().derived()); }
    };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
