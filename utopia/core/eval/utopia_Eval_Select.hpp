#ifndef UTOPIA_EVAL_SELECT_HPP
#define UTOPIA_EVAL_SELECT_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Select<InnerExpr, 1>, Traits, Backend> {
    public:
        using Expr = utopia::Select<InnerExpr, 1>;
        using Result = EXPR_TYPE(Traits, Expr);
        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&left = Eval<InnerExpr, Traits>::apply(expr.expr());

            if (result.is_alias(left)) {
                Result temp;

                left.select(expr.index(), temp);

                result.construct(std::move(temp));

            } else {
                left.select(expr.index(), result);
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Select<InnerExpr, 2>, Traits, Backend> {
    public:
        using Expr = utopia::Select<InnerExpr, 1>;
        using Result = EXPR_TYPE(Traits, InnerExpr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&mat = Eval<InnerExpr, Traits>::apply(expr.expr());

            if (result.is_alias(mat)) {
                Result temp;

                mat.select(expr.row_index(), expr.col_index(), temp);

                mat.construct(std::move(temp));

            } else {
                mat.select(expr.row_index(), expr.col_index(), result);
            }

            UTOPIA_TRACE_END(expr);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_EVAL_SELECT_HPP
