//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP
#define UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    template<class Expr, class Traits, int Backend>
    class Eval<Differentiable<Expr>, Traits, Backend> {
    public:
        EXPR_TYPE(Traits, Expr) eval(const Differentiable<Expr> &expr) {
            std::cerr << "[Warning] called eval(Differentiable<Expr>) placeholder" << std::endl;
            return Eval<Expr, Traits>::apply(expr.expr());
        }
    };

    template<class Tensor, int Order, class Traits, int Backend>
    class Eval< Differentiable<Wrapper<Tensor, Order> >, Traits, Backend> {
    public:
        inline static const Tensor &apply(const Differentiable<Wrapper<Tensor, Order> > &expr) {
            return Eval<Wrapper<Tensor, Order>, Traits>::apply(expr.expr());
        }
    };

    template<class Tensor, int Order, class Traits, int Backend>
    class Eval<Differentiable<const Wrapper<Tensor, Order> >, Traits, Backend> {
    public:
        inline static const Tensor & apply(const Differentiable<const Wrapper<Tensor, Order> > &expr)
        {
            return Eval<Wrapper<Tensor, Order>, Traits>::apply(expr.expr());
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_DIFFERENTIABLE_HPP
