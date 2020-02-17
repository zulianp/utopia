#ifndef UTOPIA_UTOPIA_EVAL_AXPY_HPP
#define UTOPIA_UTOPIA_EVAL_AXPY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    template<class LeftDerived, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Tensor<LeftDerived, Order>, Binary<Number<ScalarT>, Tensor<RightTensor, Order>, Multiplies>, Plus>, Traits, Backend> {
    public:
        using Left   = utopia::Tensor<LeftDerived, Order>;
        using Right  = utopia::Tensor<RightTensor, Order>;
        using Expr   = utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus>;
        using Result = EXPR_TYPE(Traits, Left);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&right = Eval<Right, Traits>::apply(expr.right().right());
            auto &&left  = Eval<Left, Traits>::apply(expr.left());
            auto &&alpha = expr.right().left();

            if(result.is_alias(right)) {
                result.scale(alpha);
                result.axpy(1.0, left);
            } else if(result.is_alias(left)) {
                result.axpy(alpha, right);
            } else {
                result.construct(left);
                result.axpy(
                    alpha,
                    right
                );
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    /// result = left + (right * alpha)
    template<class LeftDerived, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<
                Tensor<LeftDerived, Order>,
                Binary<
                    Tensor<RightTensor, Order>,
                    Number<ScalarT>, Multiplies>,
                    Plus
                >, Traits, Backend> {
    public:
        using Left  = utopia::Tensor<LeftDerived, Order>;
        using Right = utopia::Tensor<RightTensor, Order>;
        using Expr  = utopia::Binary<Left, Binary<Right, Number<ScalarT>, Multiplies>, Plus>;
        using Result = EXPR_TYPE(Traits, Left);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left  = Eval<Left, Traits>::apply(expr.left());
            auto &&right = Eval<Right, Traits>::apply(expr.right().left());
            auto alpha = expr.right().right();

            if(result.is_alias(right)) {
                result.scale(alpha);
                result.axpy(1.0, left);
            } else {
                if(!result.is_alias(left)) {
                    result.construct(left);
                }

                result.axpy(
                    alpha,
                    right
                );
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    /// left + alpha * identity()
    template<class LeftDerived, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Tensor<LeftDerived, 2>, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus>, Traits, Backend> {
    public:
        using Left = utopia::Tensor<LeftDerived, 2>;
        using Expr = utopia::Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus>;
        using Result = EXPR_TYPE(Traits, Left);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<Left, Traits>::apply(expr.left());
            result.shift_diag(expr.right().left().get());

            UTOPIA_TRACE_END(expr);
        }
    };

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    /// axpy specialization result = x + alpha * y
    template<class LeftDerived, int LeftOrder, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Tensor<LeftDerived, LeftOrder>,
                        Binary<Number<ScalarT>,
                             Right,
                             Multiplies
                            >,
                        Plus
                        >, Traits, Backend> {
    public:
        using Left = utopia::Tensor<LeftDerived, LeftOrder>;
        using Expr = utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&alpha = expr.right().left();
            auto &&y = Eval<Left, Traits, Backend>::apply(expr.left());
            auto &&x = Eval<Right, Traits, Backend>::apply(expr.right().right());

            if(result.is_alias(y)) {
                result.axpy(alpha, x);
            } else if(result.is_alias(x)) {
                result.scale(alpha);
                result.axpy(1.0, y);
            } else {
                result.construct(x);
                result.scale(alpha);
                result.axpy(1.0, y);
            }
            
            UTOPIA_TRACE_END(expr);
        }
    };

    /// result = l - alpha * r;
    template<class LeftDerived, int LeftOrder, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<
                    Tensor<LeftDerived, LeftOrder>,
                    Binary<
                        Number<ScalarT>,
                        Right,
                        Multiplies
                        >,
                    Minus
                    >, Traits, Backend> {
    public:
        using Left = utopia::Tensor<LeftDerived, LeftOrder>;
        using Expr = utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Minus>;
        using Result = EXPR_TYPE(Traits, Expr);
        using Scalar = typename Traits::Scalar;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&l = Eval<Left, Traits, Backend>::apply(expr.left());
            auto &&r = Eval<Right, Traits, Backend>::apply(expr.right().right());
            Scalar alpha = expr.right().left();

            if(result.is_alias(l)) {
                result.axpy(-alpha, r);
            } else if(result.is_alias(r)) {
                result.scale(-alpha);
                result.axpy(1.0, l);
            } else {
                result.construct(l);
                result.axpy(-alpha, r);
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    /// result = alpha * left + right 
    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<
                    Binary<
                        Number<ScalarT>,
                        Left,
                        Multiplies
                        >,
                    Right,
                    Plus
                    >, Traits, Backend> {
    public:
        using Expr = utopia::Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus>;
        using Result = EXPR_TYPE(Traits, Expr);
        using Scalar = typename Traits::Scalar;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            const Scalar alpha = expr.left().left();
            auto &&left  = Eval<Left, Traits, Backend>::apply(expr.left().right());
            auto &&right = Eval<Right, Traits, Backend>::apply(expr.right());

            if(result.is_alias(right)) {
                result.axpy(
                    alpha,
                    left
                );

            } else {
                if(result.is_alias(left)) {
                    result.scale(alpha);
                    result.axpy(1, right);
                } else {
                    result.assign(right);
                    result.axpy(alpha, left);
                }
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    /// result = alpha * left - right; (in order to avoid ambiguos specialization Left )
    template<typename ScalarT, class LeftDerived, int LeftOrder, class Right, class Traits, int Backend>
    class Eval<
            Binary<Binary<Number<ScalarT>, Tensor<LeftDerived, LeftOrder>, Multiplies>, Right, Minus>,
            Traits,
            Backend> {
    public:
        using Left = utopia::Tensor<LeftDerived, LeftOrder>;
        using Expr = utopia::Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Minus>;
        using Result = EXPR_TYPE(Traits, Expr);
        using Scalar = typename Traits::Scalar;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left = Eval<Left, Traits>::apply(expr.left().right());
            auto &&right = Eval<Right, Traits>::apply(expr.right());
            const ScalarT alpha = expr.left().left();

            if(result.is_alias(right)) {
                result.scale(-1.0);
                result.axpy(alpha, left);
            } else if(result.is_alias(left)) {
                result.scale(alpha);
                result.axpy(-1.0, right);
            } else {
                result.assign(left);
                result.scale(alpha);
                result.axpy(-1.0, right);
            }

            UTOPIA_TRACE_END(expr);
        }

    };

}

#endif //UTOPIA_UTOPIA_EVAL_AXPY_HPP
