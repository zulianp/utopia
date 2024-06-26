//
// Created by Patrick Zulian on 29/05/15.
//

#ifndef UTOPIA_UTOPIA_BOOLEAN_HPP
#define UTOPIA_UTOPIA_BOOLEAN_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template <class Expr>
    class Boolean : public Expression<Boolean<Expr> > {
    public:
        const Expr &expr() const { return _expr; }

        inline operator bool() const {
            // Evaluator<typename Traits<Expr>::Vector, Traits<Expr>::Backend> eval;
            // return eval.eval(expr());
            return Eval<Expr, Traits<Expr>, Traits<Expr>::Backend>::apply(expr());
        }

        Boolean(const Expr &expr) : _expr(expr) {}

    private:
        UTOPIA_STORE_CONST(Expr) _expr;
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_BOOLEAN_HPP
