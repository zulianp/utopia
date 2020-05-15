//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP
#define UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    template <class Expr, class Traits, int Backend>
    class Eval<Differentiable<Expr>, Traits, Backend> {
    public:
        EXPR_TYPE(Traits, Expr) eval(const Differentiable<Expr> &expr) {
            std::cerr << "[Warning] called eval(Differentiable<Expr>) placeholder" << std::endl;
            return Eval<Expr, Traits>::apply(expr.expr());
        }
    };

    template <class Derived, int Order, class Traits, int Backend>
    class Eval<Differentiable<Tensor<Derived, Order> >, Traits, Backend> {
    public:
        inline static const Derived &apply(const Differentiable<Tensor<Derived, Order> > &expr) {
            return Eval<Tensor<Derived, Order>, Traits>::apply(expr.expr());
        }
    };

    template <class Derived, int Order, class Traits, int Backend>
    class Eval<Differentiable<const Tensor<Derived, Order> >, Traits, Backend> {
    public:
        inline static const Derived &apply(const Differentiable<const Tensor<Derived, Order> > &expr) {
            return Eval<Tensor<Derived, Order>, Traits>::apply(expr.expr());
        }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP
