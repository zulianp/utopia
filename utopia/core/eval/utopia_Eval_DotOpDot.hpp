#ifndef UTOPIA_EVAL_DOT_OPT_DOT_HPP
#define UTOPIA_EVAL_DOT_OPT_DOT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_DotVecVecs.hpp"

namespace utopia {

    template<class X, class FunOfX, class Op>
    using DotOpDot = utopia::Binary<
                                Dot<Tensor<X, 1>, Tensor<X, 1> >,
                                Dot<FunOfX, Tensor<X, 1> >,
                                Op
                            >;

    template<class X, class FunOfX, class Op, class Traits, int Backend>
    class Eval< DotOpDot<X, FunOfX, Op>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;
        using Expr   = utopia::DotOpDot<X, FunOfX, Op>;
        using T1     = utopia::Tensor<X, 1>;

        inline static Scalar apply(const Expr &expr) {
            auto &&t1 = Eval<T1, Traits>::apply(expr.left().expr().left());
            auto &&t2 = Eval<T1, Traits>::apply(expr.left().expr().right());

            auto &&t3 = Eval<FunOfX, Traits>::apply(expr.right().expr().left());
            auto &&t4 = Eval<T1, Traits>::apply(expr.right().expr().right());

            Scalar left, right;
            EvalDots<X, Backend>::apply(
                t1, t2, left,
                t3, t4, right
            );

            return expr.operation().apply(left, right);
        }
    };

}

#endif //UTOPIA_EVAL_DOT_OPT_DOT_HPP
