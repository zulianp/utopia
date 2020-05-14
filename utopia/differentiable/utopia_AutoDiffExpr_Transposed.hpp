#ifndef UTOPIA_AUTO_DIFF_EXPR_TRANSPOSED_HPP
#define UTOPIA_AUTO_DIFF_EXPR_TRANSPOSED_HPP

#include "utopia_AutoDiffExpr.hpp"

namespace utopia {

    template<class Expr>
    class AutoDiffExpr<Transposed<Expr>, 1> {
    public:
        using Type = utopia::Transposed<typename AutoDiffExpr<Expr>::Type>;

        inline static UTOPIA_STORE_CONST(Type) make(const Transposed<Expr> &expr)
        {
            return transpose(AutoDiffExpr<Expr>::make(expr.expr()));
        }
    };


    template<int Order>
    class Simplify< Transposed< Factory<Identity, Order> > > {
    public:
        typedef utopia::Factory<Identity, Order> Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Transposed< Factory<Identity, Order> > &expr)
        {
            return Type(size(expr));
        }
    };

    template<int Order>
    class Simplify< Transposed< Factory<Zeros, Order> > > {
    public:
        typedef utopia::Factory<Zeros, Order> Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Transposed< Factory<Zeros, Order> > &expr)
        {
            return Type(size(expr));
        }
    };

    template<class Expr>
    class Simplify< Transposed< Transposed<Expr> > > {
    public:
        using Type = typename Simplify<Expr>::Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Transposed< Transposed<Expr> > &expr)
        {
            return Simplify<Expr>::make(expr.expr().expr());
        }
    };
}

#endif //UTOPIA_AUTO_DIFF_EXPR_TRANSPOSED_HPP
