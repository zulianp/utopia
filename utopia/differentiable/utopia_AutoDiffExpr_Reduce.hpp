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

        using DiffLeft = utopia::AutoDiffExpr<Left>;
        using DiffRight = utopia::AutoDiffExpr<Right>;

        using DLeft = typename DiffLeft::Type;
        using DRight = typename DiffRight::Type;

        typedef utopia::Binary<Multiply<Transposed<DLeft>,  Right>,
                               Multiply<Transposed<DRight>, Left>,
                               Plus> ComplexType;

        using Sim = utopia::Simplify<ComplexType>;
        using Type = typename Sim::Type;

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
