#ifndef UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP
#define UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP


#include "utopia_Expressions.hpp"
#include "utopia_AutoDiffExpr.hpp"
#include "utopia_AutoDiffExpr_Diag.hpp"

namespace utopia {

    template<class Inner>
    class AutoDiffExpr< Trace<Inner>, 1>  {
    public:
        typedef utopia::Trace<Inner> Expr;

        typedef utopia::AutoDiffExpr< Diag<Inner> > Diff;
        typedef utopia::Diag<typename Diff::Type> ComplexType;
        typedef utopia::Simplify<ComplexType> Sim;

        typedef typename Sim::Type Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Expr &expr)
        {
            // unused
            //const auto &e = expr.expr();
            return Sim::make( diag( Diff::make( diag(expr.expr())) ) );
        }
    };
}

#endif //UTOPIA_AUTO_DIFF_EXPR_TRACE_HPP