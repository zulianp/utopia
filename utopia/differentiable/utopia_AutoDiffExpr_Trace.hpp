#ifndef UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP
#define UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP

#include "utopia_AutoDiffExpr.hpp"
#include "utopia_AutoDiffExpr_Diag.hpp"
#include "utopia_Expressions.hpp"

namespace utopia {

    template <class Inner>
    class AutoDiffExpr<Trace<Inner>, 1> {
    public:
        using Expr = utopia::Trace<Inner>;

        using Diff = utopia::AutoDiffExpr<Diag<Inner> >;
        using ComplexType = utopia::Diag<typename Diff::Type>;
        using Sim = utopia::Simplify<ComplexType>;

        using Type = typename Sim::Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Expr &expr) {
            // unused
            // const auto &e = expr.expr();
            return Sim::make(diag(Diff::make(diag(expr.expr()))));
        }
    };
}  // namespace utopia

#endif  // UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP