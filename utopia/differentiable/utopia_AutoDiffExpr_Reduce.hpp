#ifndef UTOPIA_AUTO_DIFF_REDUCE_HPP
#define UTOPIA_AUTO_DIFF_REDUCE_HPP

#include "utopia_Expressions.hpp"
#include "utopia_AutoDiffExpr.hpp"


namespace utopia {

    template<class Left, class Right>
    class AutoDiffExpr< Reduce<
                             Binary<Left, Right, EMultiplies>,
                             Plus>, 1>  {
    public:
        typedef utopia::Reduce< Binary<Left, Right, EMultiplies>, Plus> Expr;

        typedef utopia::AutoDiffExpr<Left>  DiffLeft;
        typedef utopia::AutoDiffExpr<Right> DiffRight;

        typedef typename DiffLeft::Type DLeft;
        typedef typename DiffRight::Type DRight;

        typedef utopia::Binary<Multiply<Transposed<DLeft>,  Right>,
                               Multiply<Transposed<DRight>, Left>,
                               Plus> ComplexType;

        typedef utopia::Simplify<ComplexType> Sim;
        typedef typename Sim::Type Type;

        inline static UTOPIA_STORE_CONST(Type) make(const Expr &expr)
        {
            const auto &e = expr.expr();
            return Sim::make(
                        transpose(DiffLeft::make(e.left())) * e.right() +
                        transpose(DiffRight::make(e.right())) * e.left());
        }
    };
}

#endif //UTOPIA_AUTO_DIFF_REDUCE_HPP
