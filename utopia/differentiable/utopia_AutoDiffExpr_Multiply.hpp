#ifndef UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP
#define UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP

#include "utopia_AutoDiffExpr.hpp"
#include "utopia_Differentiable.hpp"
#include "utopia_Simplify.hpp"

namespace utopia {

    template <class Left, class Right>
    class AutoDiffExpr<Multiply<Left, Right>, 1> {
    public:
        using DiffLeft = typename utopia::AutoDiffExpr<Left>;
        using DiffRight = typename utopia::AutoDiffExpr<Right>;

        using DLeft = typename DiffLeft::Type;
        using DRight = typename DiffRight::Type;

        typedef utopia::Binary<utopia::Multiply<DLeft, Right>, utopia::Multiply<Left, DRight>, Plus> ComplexType;

        using Type = typename utopia::Simplify<ComplexType>::Type;

        static UTOPIA_STORE_CONST(Type) make(const Multiply<Left, Right> &expr) {
            return utopia::Simplify<ComplexType>::make(DiffLeft::make(expr.left()) * expr.right() +
                                                       expr.left() * DiffRight::make(expr.right()));
        }
    };
}  // namespace utopia

#endif  // UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP
