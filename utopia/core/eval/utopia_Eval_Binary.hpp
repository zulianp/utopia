//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Operators.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template<class Result>
    class EvalBinaryAux {
    public:
        using Scalar = typename Traits<Result>::Scalar;

        template<class Left, class Right, int Order>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Plus &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.axpy(1.0, right.derived());
        }

        template<class Left, class Right, int Order>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Minus &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.axpy(-1.0, right.derived());
        }

        template<class Left, class Right>
        static void apply(const Tensor<Left, 2> &left, const Tensor<Right, 1> &right, const Multiplies &, Result &result)
        {
            left.derived().multiply(right.derived(), result);
        }

        template<class Left, class Right>
        static void apply(const Tensor<Left, 2> &left, const Tensor<Right, 2> &right, const Multiplies &, Result &result)
        {
            left.derived().multiply(right.derived(), result);
        }

        template<class Left, class Right, int Order>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Min &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_min(right.derived());
        }

        template<class Left, class Right, int Order>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Max &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_max(right.derived());
        }
 
        template<class Left, class Right>
        static void apply(Left &&left, const Number<Right> &right, const Plus &, Result &result)
        {
            // result.construct(std::forward<Left>(left));
            // result.axpy(1.0, right);
            assert(false);
        }

        template<class Left, class Right>
        static void apply(Left &&left, const Number<Right> &right, const Minus &, Result &result)
        {
            // result.construct(std::forward<Left>(left));
            // result.axpy(-1.0, right);
            assert(false);
        }

        template<class Left, class Right>
        static void apply(Left &&left, const Number<Right> &right, const Multiplies &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.scale(right);
        }

        template<class Left>
        static void apply(Left &&left, const Scalar &right, const Multiplies &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.scale(right);
        }

        template<class Left, class Right>
        static void apply(const Number<Left> &left, Right &&right, const Multiplies &, Result &result)
        {
            result.construct(std::forward<Right>(right));
            result.scale(left);
        }

        template<class Left, class Right>
        static void apply(Left &&left, const Number<Right> &right, const Min &, Result &result)
        {
            // result.construct(std::forward<Left>(left));
            // result.e_min(right);
            assert(false);
        }

        template<class Left, class Right>
        static void apply(Left &&left, const Number<Right> &right, const Max &, Result &result)
        {
            // result.construct(std::forward<Left>(left));
            // result.e_max(right);

            assert(false);
        }
    };

    template<class ScalarT, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Number<ScalarT>, Right, Operation>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Binary<Number<ScalarT>, Right, Operation> >::Type Result;
        inline static Result apply(const Binary<Number<ScalarT>, Right, Operation> &expr)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            EvalBinaryAux<Result>::apply(
                Eval<Right, Traits>::apply(expr.right()),
                expr.left(),
                expr.operation(),
                result
            );


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class ScalarT, class Left, class Operation, class Traits, int Backend>
    class Eval<Binary<Left, Number<ScalarT>, Operation>, Traits, Backend> {
    public:
        using Expr = utopia::Binary<Left, Number<ScalarT>, Operation>;

        typedef typename TypeAndFill<Traits, Expr>::Type Result;
        inline static Result apply(const Expr &expr)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            EvalBinaryAux<Result>::apply(
                Eval<Left, Traits>::apply(expr.left()),
                expr.right(),
                expr.operation(),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Number<Left>, Number<Right>, Operation>, Traits, Backend> {
    public:

        inline static auto apply(const Binary<Number<Left>, Number<Right>, Operation> &expr) -> decltype(Left() + Right())
        {
            Left l = expr.left();
            Right r = expr.right();
            return expr.operation().apply(l, r);
        }
    };

 

    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Left, Right, Operation>, Traits, Backend> {
    public:
        typedef typename utopia::TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type Result;

        inline static Result apply(const Binary<Left, Right, Operation> &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            EvalBinaryAux<Result>::apply(
                Eval<Left,  Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right()),
                expr.operation(),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class ScalarT, class Traits, int Backend>
    class Eval<Binary<Left, Number<ScalarT>, Multiplies>, Traits, Backend> {
    public:

        typedef typename TypeAndFill<Traits, Binary<Left, Number<ScalarT>, Multiplies> >::Type Result;
        inline static Result apply(const Binary<Left, Number<ScalarT>, Multiplies> &expr)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct(
                Eval<Left, Traits>::apply(expr.left())
            );

            result.scale(expr.right());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };


    template<class Left, class Right, class Traits, int Backend>
    class Eval<OuterProduct<Left, Right>, Traits, Backend> {
    public:
        typedef utopia::OuterProduct<Left, Right> Expr;

        inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr) {
            EXPR_TYPE(Traits, Expr) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).kronecker_product(
                result,
                Eval<Left, Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right())
                );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
