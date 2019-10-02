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
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &&t1 = Eval<T1, Traits>::apply(expr.left().expr().left());
            auto &&t2 = Eval<T1, Traits>::apply(expr.left().expr().right());

            auto &&t3 = Eval<FunOfX, Traits>::apply(expr.right().expr().left());
            auto &&t4 = Eval<T1, Traits>::apply(expr.right().expr().right());

            Scalar left, right;
            EvalDots<X, Backend>::apply(
                t1, t2, left,
                t3, t4, right
            );

            const Scalar ret = expr.operation().apply(left, right);

            UTOPIA_TRACE_END_SPECIALIZED(expr);
            return ret;
        }
    };

    template<class X, typename T, class FunOfX, class Op, class LeftOp, class RightOp>
    using AlphaDotOpBetaDot = utopia::Binary<
                                Binary<Number<T>, Dot<Tensor<X, 1>, Tensor<X, 1> >, LeftOp>,
                                Binary<Number<T>, Dot<FunOfX, Tensor<X, 1> >, RightOp>,
                                Op
                            >;

    template<class X, typename T, class FunOfX, class Op, class LeftOp, class RightOp, class Traits, int Backend>
    class Eval< AlphaDotOpBetaDot<X, T, FunOfX, Op, LeftOp, RightOp>, Traits, Backend> {
    public:
        using Expr = utopia::AlphaDotOpBetaDot<X, T, FunOfX, Op, LeftOp, RightOp>;

        using Scalar = typename Traits::Scalar;
        using T1     = utopia::Tensor<X, 1>;

        inline static Scalar apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &&t1 = Eval<T1, Traits>::apply(expr.left().right().expr().left());
            auto &&t2 = Eval<T1, Traits>::apply(expr.left().right().expr().right());

            auto &&t3 = Eval<FunOfX, Traits>::apply(expr.right().right().expr().left());
            auto &&t4 = Eval<T1, Traits>::apply(expr.right().right().expr().right());

            Scalar alpha = expr.left().left();
            Scalar beta  = expr.right().left();

            auto &&left_op = expr.left().operation();
            auto &&right_op = expr.right().operation();

            Scalar left, right;
            EvalDots<X, Backend>::apply(
                t1, t2, left,
                t3, t4, right
            );

            const Scalar ret = expr.operation().apply(
                    left_op.apply(alpha, left), 
                    right_op.apply(beta, right)
                );

            UTOPIA_TRACE_END_SPECIALIZED(expr);
            return ret;
        }
    };
}

#endif //UTOPIA_EVAL_DOT_OPT_DOT_HPP
