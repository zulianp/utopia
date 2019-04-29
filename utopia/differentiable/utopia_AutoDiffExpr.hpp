#ifndef UTOPIA_AUTO_DIFF_EXPR_HPP
#define UTOPIA_AUTO_DIFF_EXPR_HPP

#include "utopia_Differentiable.hpp"
#include "utopia_Simplify.hpp"

namespace utopia {
    template<class Expr, int IsDifferentiable = is_differentiable<Expr>() >
    class AutoDiffExpr {};

    template<class Expr>
    class AutoDiffExpr<Differentiable<Expr>, 1> {
    public:
        typedef utopia::Factory<Identity, Expr::Order * 2> Type;

        static Type make(const Differentiable<Expr> &expr)
        {
            auto s = size(expr.expr());
            return Type(
                        kron_prod(s, s)
                    );
        }
    };

    template<class Expr>
    class AutoDiffExpr<Expr, 0> {
    public:
        static const int Order = Expr::Order;

        typedef UTOPIA_SCALAR(Expr) Scalar;
        typedef utopia::Factory<Zeros, Order> Type;

        static Type make(const Differentiable<Expr> &expr)
        {
            return Type(size(expr));
        }
    };

    template<class InnerExpr, class Operation>
    class AutoDiffUnary {};

    template<class InnerExpr>
    class AutoDiffUnary<InnerExpr, Pow2> {
    public:
        typedef typename InnerExpr::Scalar Scalar;
        //[new backend concept]
        typedef Diag< utopia::Binary< Number<Scalar>, InnerExpr, utopia::Multiplies> > Type;
        // typedef Unary< utopia::Binary< Number<Scalar>, InnerExpr, utopia::Multiplies>, DiagOp> Type;

        static UTOPIA_STORE_CONST(Type) make(const InnerExpr &expr, const Pow2 &)
        {
            return diag(Scalar(2.0) * expr);
        }
    };

    //d x^(1/2)  = 0.5 * x^(-1/2)
    template<class InnerExpr>
    class AutoDiffUnary<InnerExpr, Sqrt> {
    public:
        typedef typename InnerExpr::Scalar Scalar;
        typedef utopia::Unary< Unary<InnerExpr, Sqrt>, Reciprocal<Scalar> > Type;

        static UTOPIA_STORE_CONST(Type) make(const InnerExpr &expr, const Sqrt &)
        {
            return Scalar(.5) / sqrt(expr);
        }
    };

    template<class InnerExpr, class Operation>
    class AutoDiffExpr< Unary<InnerExpr, Operation>, 1> {
    public:
        typedef utopia::AutoDiffUnary<InnerExpr, Operation> Diff;
        typedef typename Diff::Type Left;
        typedef typename AutoDiffExpr<InnerExpr>::Type Right;

        typedef utopia::Multiply<Left, Right> ComplexType;
        // Type;
        typedef typename utopia::Simplify<ComplexType>::Type Type;


        static UTOPIA_STORE_CONST(Type) make(const Unary<InnerExpr, Operation> &expr)
        {
            return 	utopia::Simplify<ComplexType>::make( ComplexType(

                                                                     Diff::make(expr.expr(), expr.operation()),
                                                                     AutoDiffExpr<InnerExpr>::make(
                                                                                                   expr.expr()
                                                                                                   )
                                                                     ));
        }
    };

    template<class Derived>
    inline typename AutoDiffExpr<Derived>::Type derivative(const Expression<Derived> &expr)
    {
        static_assert(is_differentiable<Derived>(), "must be differentiable: use independent_variable");
        return AutoDiffExpr<Derived>::make(expr.derived());
    }
}


#endif //UTOPIA_AUTO_DIFF_EXPR_HPP
