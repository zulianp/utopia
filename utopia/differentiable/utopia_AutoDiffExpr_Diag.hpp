#ifndef UTOPIA_AUTO_DIFF_EXPR_DIAG_HPP
#define UTOPIA_AUTO_DIFF_EXPR_DIAG_HPP

#include "utopia_Expressions.hpp"
#include "utopia_AutoDiffExpr.hpp"

namespace utopia {

    template<class Inner, int Order = Inner::Order>
    class AutoDiffExprDiag {};


    template<class Inner>
    class AutoDiffExprDiag<Inner, 2> {
    public:
        //[new backend concept]
        typedef utopia::Diag<typename AutoDiffExpr<Inner>::Type> Type;
        // typedef utopia::Unary<typename AutoDiffExpr<Inner>::Type, DiagOp> Type;


        inline static UTOPIA_STORE_CONST(Type) make(const Inner &expr)
        {
            return diag(AutoDiffExpr<Inner>::make(expr));
        }
    };

    template<class Inner>
    class AutoDiffExprDiag<Inner, 1> {
    public:
        static_assert(Inner::Order > 1, "not supported");
        typedef Inner Type;
    };

    template<class Inner>
    class AutoDiffExpr< Diag<Inner>, 1>  {
    public:
        //[new backend concept]
        typedef utopia::Diag<Inner> Expr;
        // typedef utopia::Unary<Inner, DiagOp> Expr;


        typedef typename AutoDiffExprDiag<Inner>::Type Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Expr &expr)
        {
            const auto &e = expr.expr();
            return AutoDiffExprDiag<Inner>::make(expr.expr());

        }
    };

}

#endif //UTOPIA_AUTO_DIFF_EXPR_DIAG_HPP
