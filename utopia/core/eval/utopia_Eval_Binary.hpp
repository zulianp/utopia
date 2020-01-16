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
        static void apply(Left &&left, Right &&right, const Plus &, Result &result)
        {
            if(result.is_alias(left)) {
                result.axpy(1.0, right);
            } else if(result.is_alias(right)) {
                result.axpy(1.0, left);
            } else {
                if(std::is_rvalue_reference<Right>::value) {
                    result.construct(std::forward<Right>(right));
                    result.axpy(1.0, left);
                } else {
                    result.construct(std::forward<Left>(left));
                    result.axpy(1.0, right);
                }
            }
        }

        ///////////////////////////// MINUS /////////////////////////////

        template<class Left, class Right>
        static void apply(Left &&left, Right &&right, const Minus &, Result &result)
        {
            if(result.is_alias(left)) {
                result.axpy(-1.0, right);
                return;
            }

            if(result.is_alias(right)) {
                //necessary for petsc
                result.axpy(-1.0, left);
                result.scale(-1.0);
                return;
            }

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
        using Expr = utopia::Binary<Left, Right, Operation>;
        using Result = EXPR_TYPE(Traits, Expr);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            apply_aux(
                Eval<Left,  Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right()),
                expr.operation(),
                result
            );

            UTOPIA_TRACE_END(expr);
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

            apply(
                Eval<TensorT, Traits>::apply(expr.right().left()),
                Eval<TensorT, Traits>::apply(expr.right().right()),
                op,
                result
            );

            UTOPIA_TRACE_END(expr);
        }

        template<class TL, class TR, class Res>
        inline static void apply(TL &&l, TR &&r, const Op &op, Res &result)
        {
            if(r.is_alias(result)) {
                if(!utopia::is_commutative<Op>::value) {
                    if(std::is_same<Divides, Op>::value) {
                        result.reciprocal(1.0);
                        result.e_mul(l);
                    } else {
                        assert(false && "IMPLEMENT ME");
                        std::cerr << "[Error] EvalAssignElWise: case not implemented " << std::endl;
                        Utopia::Abort();
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

    template<typename T, class Op, class Traits, int Backend>
    class Eval< Binary<Tensor<T, 1>, Tensor<T, 1>, Op>, Traits, Backend> {
    public:
        using Expr = utopia::Binary<Tensor<T, 1>, Tensor<T, 1>, Op>;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, T)

        inline static void apply(const Expr &expr, T &result)
        {
            EvalAssignElWise<Tensor<T, 1>, Op, Traits>::apply(
                Eval<Tensor<T, 1>, Traits>::apply(expr.left()),
                Eval<Tensor<T, 1>, Traits>::apply(expr.right()),
                expr.operation(),
                result.derived()
            );
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

    template<class T, class Traits, int Backend>
    class Eval<TensorAssignElWise<T, Divides>, Traits, Backend> {
    public:
        
        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            EvalAssignElWise<Tensor<T, 1>, Divides, Traits>::apply(expr);
            UTOPIA_TRACE_END(expr);
        }
    };

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
        using Result = EXPR_TYPE(Traits, Expr);

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left  = Eval<Left, Traits>::apply(expr.left());
            auto &&right = Eval<Right, Traits>::apply(expr.right());

            EvalKroneckerProduct<Result, typename Traits::Vector>::apply(left, right, result);

            UTOPIA_TRACE_END(expr);
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_BINARY_HPP_HPP
