#ifndef UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Operators.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Eval_KroneckerProduct.hpp"
#include "utopia_Tracer.hpp"

namespace utopia {

    template<class Result>
    class EvalBinaryAux {
    public:
        template<class Left, class Right, class Operation, class Res, int Order>
        static void apply(Left &&l, Right &&r, const Operation &op, Tensor<Res, Order> &result)
        {
            EvalBinaryAux<Tensor<Res, Order>>::apply(
                std::forward<Left>(l),
                std::forward<Right>(r),
                op,
                result.derived()
            );
        }

        template<class Left, class Right, class Operation, class Res>
        static void apply(Left &&l, Right &&r, const Operation &op, Number<Res> &result)
        {
            EvalBinaryAux<Number<Res>>::apply(
                std::forward<Left>(l),
                std::forward<Right>(r),
                op,
                result
            );
        }
    };

    template<class Result, int Order>
    class EvalBinaryAux<Tensor<Result, Order>> {
    public:
        using Scalar = typename Traits<Result>::Scalar;

        ///////////////////////////// PLUS /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Plus &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.axpy(1.0, right.derived());
        }

        ///////////////////////////// MINUS /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Minus &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.axpy(-1.0, right.derived());
        }

        ///////////////////////////// EMULTIPLIES /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const EMultiplies &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_mul(right.derived());
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

        ///////////////////////////// EDIVIDES /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Divides &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_div(right.derived());
        }

        template<class Left, class Right>
        static void apply(Left &&left, const Scalar &right, const Divides &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_div(right);
        }

        ///////////////////////////// MULTIPLIES /////////////////////////////

        template<class Left>
        static void apply(Left &&left, const Scalar &right, const Multiplies &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.scale(right);
        }

        template<class Right>
        static void apply(const Scalar &left, Right &&right, const Multiplies &, Result &result)
        {
            result.construct(std::forward<Right>(right));
            result.scale(left);
        }

        ///////////////////////////// MIN /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Min &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_min(right.derived());
        }

        template<class Left>
        static void apply(Left &&left, const Scalar &right, const Min &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_min(right);
        }

        template<class Right>
        static void apply(const Scalar &l, Right &&right, const Min &, Result &result)
        {
            result.construct(std::forward<Right>(right));
            result.e_min(l);
        }

        ///////////////////////////// MAX /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, const Tensor<Right, Order> &right, const Max &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_max(right.derived());
        }


        template<class Right>
        static void apply(const Scalar &l, Right &&right, const Max &, Result &result)
        {
            result.construct(std::forward<Right>(right));
            result.e_max(l);
        }

        template<class Left>
        static void apply(Left &&left, const Scalar &r, const Max &, Result &result)
        {
            result.construct(std::forward<Left>(left));
            result.e_max(r);
        }
    };

    template<>
    class EvalBinaryAux<double> {
    public:
        template<class Op>
        static void apply(const double &left, const double &right, const Op &op, double &result)
        {
             result = op.template apply<double>(left, right);
        }
    };

    template<typename Scalar>
    class EvalBinaryAux<Number<Scalar>> {
    public:
        template<class Op>
        static void apply(const Number<Scalar> &left, const Number<Scalar> &right, const Op &op, Number<Scalar> &result)
        {
             result.set(op.template apply<double>(left.get(), right.get()));
        }
    };

    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<Binary<Left, Right, Operation>, Traits, Backend> {
    public:
        typedef typename utopia::TypeAndFill<Traits, Binary<Left, Right, Operation> >::Type Result;

        inline static Result apply(const Binary<Left, Right, Operation> &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            apply_aux(
                Eval<Left,  Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right()),
                expr.operation(),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }

        template<class L, class R, class Res, int Order>
        inline static void apply_aux(L &&l, R &&r, const Operation &op, Tensor<Res, Order> &result)
        {
            EvalBinaryAux<Tensor<Res, Order>>::apply(
                std::forward<L>(l),
                std::forward<R>(r),
                op,
                result.derived()
            );
        }

        template<class L, class R, class Res>
        inline static void apply_aux(L &&l, R &&r, const Operation &op, Number<Res> &result)
        {
            EvalBinaryAux<Number<Res>>::apply(
                std::forward<L>(l),
                std::forward<R>(r),
                op,
                result
            );
        }
    };

    template<typename T, class Op>
    using TensorAssignElWise =  Assign<
                Tensor<T, 1>,
                Binary<Tensor<T, 1>, Tensor<T, 1>, Op>>;

    template<class TensorT, class Op, class Traits>
    class EvalAssignElWise {
    public:
        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &result = Eval<TensorT, Traits>::apply(expr.left());

            auto &&op = expr.right().operation();

            apply_aux(
                Eval<TensorT, Traits>::apply(expr.right().left()),
                Eval<TensorT, Traits>::apply(expr.right().right()),
                op,
                result
            );

            UTOPIA_TRACE_END(expr);
        }

        template<class TL, class TR, class Res>
        inline static void apply_aux(TL &&l, TR &&r, const Op &op, Res &result)
        {
            if(r.same_object(result)) {
                if(!utopia::is_commutative<Op>::value) {
                    if(std::is_same<Divides, Op>::value) {
                        result.reciprocal(1.0);
                        result.e_mul(r);
                    } else {
                        assert(false && "IMPLEMENT ME");
                    }

                } else {
                    EvalBinaryAux<TensorT>::apply(
                        std::forward<TR>(r),
                        std::forward<TL>(l),
                        op,
                        result.derived()
                    );
                }
            } else {
                EvalBinaryAux<TensorT>::apply(
                    std::forward<TL>(l),
                    std::forward<TR>(r),
                    op,
                    result.derived()
                );
            }
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<TensorAssignElWise<T, EMultiplies>, Traits, Backend> {
    public:

        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            EvalAssignElWise<Tensor<T, 1>, EMultiplies, Traits>::apply(expr);
            UTOPIA_TRACE_END(expr);
        }
    };

    // template<class T, class Traits, int Backend>
    // class Eval<TensorAssignElWise<T, Divides>, Traits, Backend> {
    // public:
        
    //     template<class Expr>
    //     inline static void apply(const Expr &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);
    //         EvalAssignElWise<Tensor<T, 1>, Divides, Traits>::apply(expr);
    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    template<class T, class Traits, int Backend>
    class Eval<TensorAssignElWise<T, Min>, Traits, Backend> {
    public:
        
        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            EvalAssignElWise<Tensor<T, 1>, Min, Traits>::apply(expr);
            UTOPIA_TRACE_END(expr);
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<TensorAssignElWise<T, Max>, Traits, Backend> {
    public:
        
        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            EvalAssignElWise<Tensor<T, 1>, Max, Traits>::apply(expr);
            UTOPIA_TRACE_END(expr);
        }
    };
   
    template<class Left, class Right, class Traits, int Backend>
    class Eval<OuterProduct<Left, Right>, Traits, Backend> {
    public:
        typedef utopia::OuterProduct<Left, Right> Expr;

        inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr) {
            EXPR_TYPE(Traits, Expr) result;

            UTOPIA_TRACE_BEGIN(expr);

            auto left  = Eval<Left, Traits>::apply(expr.left());
            auto right = Eval<Left, Traits>::apply(expr.right());

            EvalKroneckerProduct<decltype(result), decltype(left)>::apply(left, right, result);

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
