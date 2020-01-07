#ifndef UTOPIA_UTOPIA_EVAL_TERNARY_HPP
#define UTOPIA_UTOPIA_EVAL_TERNARY_HPP

#include "utopia_Eval_Binary.hpp"


namespace utopia {

    template<class T1, class T2, class T3, int Order, class Op1, class Op2>
    using TernaryExpr = utopia::Binary<utopia::Binary<Tensor<T1, Order>, Tensor<T2, Order>, Op1>, Tensor<T3, Order>, Op2>;

    template<class T1, class T2, class T3, int Order, class Op1, class Op2, class Traits, int Backend>
    class Eval<TernaryExpr<T1, T2, T3, Order, Op1, Op2>, Traits, Backend> {
    public:
        using Expr   = utopia::TernaryExpr<T1, T2, T3, Order, Op1, Op2>;
        using Result = T1;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&t1 = Eval<Tensor<T1, Order>, Traits>::apply(expr.left().left());
            auto &&t2 = Eval<Tensor<T2, Order>, Traits>::apply(expr.left().right());
            auto &&t3 = Eval<Tensor<T3, Order>, Traits>::apply(expr.right());

            const auto op1 = expr.left().operation();
            const auto op2 = expr.operation();

            const static bool order_matters = !eval_order_changable<Op1, Op2>::value;

            if(result.is_alias(t3)) {
                if(!order_matters) {
                    apply_aux(result, t3, t2, t1, op1, op2);
                } else {
                    //Temporary is created here
                    auto temp = t3;
                    apply_aux(result, t1, t2, temp, op1, op2);
                }

            } else {
                apply_aux(result, t1, t2, t3, op1, op2);
            }

            UTOPIA_TRACE_END(expr);
        }

    private:
        template<class AuxR, class AuxT1, class AuxT2, class AuxT3, class AuxOp1, class AuxOp2>
        inline static void apply_aux(AuxR &result, AuxT1 &&t1, AuxT2 &&t2, AuxT3 &&t3, const AuxOp1 &op1, const AuxOp2 &op2)
        {
            EvalBinaryAux<AuxR>::apply(t1, t2, op1, result);
            EvalBinaryAux<AuxR>::apply(result, t3, op2, result);
        }
    };

    template<class T1, class T2, class T3, int Order, class Op1, class Op2, class Traits, int Backend>
    using EvalTernary = utopia::Eval<TernaryExpr<T1, T2, T3, Order, Op1, Op2>, Traits, Backend>;

    template<class Result, class T1, class T2, class T3, int Order, class Op1, class Op2, class Traits, int Backend>
    class Eval<Assign<Tensor<Result, Order>, TernaryExpr<T1, T2, T3, Order, Op1, Op2>>, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Tensor<Result, Order>, TernaryExpr<T1, T2, T3, Order, Op1, Op2>>;

        inline static void apply(const Expr &expr) {
            EvalTernary<T1, T2, T3, Order, Op1, Op2, Traits, Backend>::apply(expr.right(), expr.left().derived());
        }
    };


    //result = a op1 (alpha * b) op2 (beta * c);
    //Assign<Vec, Minus<Plus<Vec, Multiplies<Number, Vec>>, Multiplies<Number, Vec>>>
    template<class V, typename T, class Op1, class Op2>
    using ScaledTernaryExpr = utopia::
                                Binary<
                                    Binary<
                                        Tensor<V, 1>,
                                        Binary<
                                            Number<T>,
                                            Tensor<V, 1>,
                                            Multiplies
                                            >,
                                        Op1
                                        >,
                                    Binary<
                                        Number<T>,
                                        Tensor<V, 1>,
                                        Multiplies
                                        >,
                                    Op2>;

    template<class V, typename T, class Traits, int Backend>
    class Eval< ScaledTernaryExpr<V, T, Plus, Minus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<V, 1>;
        using Scalar = typename Traits::Scalar;
        using Expr   = utopia::ScaledTernaryExpr<V, T, Plus, Minus>;
        using Result = V;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);
            //result = (a + (alpha * b)) - (beta * c);
            auto &a = Eval<T1, Traits>::apply(expr.left().left());
            const Scalar alpha = expr.left().right().left();
            auto &b = Eval<T1, Traits>::apply(expr.left().right().right());
            const Scalar beta = expr.right().left();
            auto &c = Eval<T1, Traits>::apply(expr.right().right());

            if(result.is_alias(b)) {
                result.scale(alpha);
                result.axpy(1.0, a);
                result.axpy(-beta, c);
            } else if(result.is_alias(c)) {
                result.scale(-beta);
                result.axpy(alpha, b);
                result.axpy(1.0, a);
            } else if(result.is_alias(a)) {
                result.axpy(alpha, b);
                result.axpy(-beta, c);
            } else {
                result = a;
                result.axpy(alpha, b);
                result.axpy(-beta, c);
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    template<class V, typename T, class Traits, int Backend>
    class Eval< Assign<Tensor<V, 1>, ScaledTernaryExpr<V, T, Plus, Minus>>, Traits, Backend> {
    public:
        static void apply(const Assign<Tensor<V, 1>, ScaledTernaryExpr<V, T, Plus, Minus>> &expr)
        {
            Eval< ScaledTernaryExpr<V, T, Plus, Minus> >::apply(expr.right(), expr.left().derived());
        }
    };


    template<class V, typename T, class Traits, int Backend>
    class Eval< ScaledTernaryExpr<V, T, Minus, Plus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<V, 1>;
        using Scalar = typename Traits::Scalar;
        using Expr = utopia::ScaledTernaryExpr<V, T, Minus, Plus>;
        using Result = V;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);
            //result = (a - (alpha * b)) + (beta * c);
            auto &a = Eval<T1, Traits>::apply(expr.left().left());
            const Scalar alpha = expr.left().right().left();
            auto &b = Eval<T1, Traits>::apply(expr.left().right().right());
            const Scalar beta = expr.right().left();
            auto &c = Eval<T1, Traits>::apply(expr.right().right());

            if(result.is_alias(b)) {
                result.scale(-alpha);
                result.axpy(1.0, a);
                result.axpy(beta, c);
            } else if(result.is_alias(c)) {
                result.scale(beta);
                result.axpy(-alpha, b);
                result.axpy(1.0, a);
            } else if(result.is_alias(a)) {
                result.axpy(-alpha, b);
                result.axpy(beta, c);
            } else {
                result = a;
                result.axpy(-alpha, b);
                result.axpy(beta, c);
            }

            UTOPIA_TRACE_END(expr);
        }
    };


    template<class V, typename T, class Traits, int Backend>
    class Eval< Assign<Tensor<V, 1>, ScaledTernaryExpr<V, T, Minus, Plus>>, Traits, Backend> {
    public:
        static void apply(const Assign<Tensor<V, 1>, ScaledTernaryExpr<V, T, Minus, Plus>> &expr)
        {
            Eval< ScaledTernaryExpr<V, T, Minus, Plus> >::apply(expr.right(), expr.left().derived());
        }
    };

    //Assign<PetscMatrix, Plus<PetscMatrix, Multiply<PetscMatrix, PetscMatrix>>>
    template<class M, class Op, class Traits, int Backend>
    class Eval< Binary<M, Multiply<M, M>, Op>, Traits, Backend> {
    public:
        using Expr   = utopia::Binary<M, Multiply<M, M>, Op>;
        using Result = EXPR_TYPE(Traits, M);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);
            // result = A op (B * C)
            auto &&A = Eval<M, Traits>::apply(expr.left());
            auto &&B = Eval<M, Traits>::apply(expr.right().left());
            auto &&C = Eval<M, Traits>::apply(expr.right().right());
            auto && op = expr.operation();

            if(result.is_alias(A) || result.is_alias(B) || result.is_alias(C)) {
                typename Unwrap<M>::Type temp;
                B.multiply(C, temp);
                apply_aux(A, std::move(temp), op, result);
            } else {
                B.multiply(C, result);
                apply_aux(A, result, op, result);
            }

            UTOPIA_TRACE_END(expr);
        }

        template<class MLeft, class MRight, class MResult>
        inline static void apply_aux(const MLeft &left, MRight &&right, const Op &op, MResult &result)
        {
            EvalBinaryAux<MResult>::apply(left, std::forward<MRight>(right), op, result);
        }

    };


    //Assign<PetscMatrix, Plus<PetscMatrix, Multiply<PetscMatrix, PetscMatrix>>>
    template<class M, class Op, class Traits, int Backend>
    class Eval< Assign<M, Binary<M, Multiply<M, M>, Op>>, Traits, Backend> {
    public:
        static void apply(const Assign<M, Binary<M, Multiply<M, M>, Op>> &expr)
        {
            Eval<Binary<M, Multiply<M, M>, Op>, Traits, Backend>::apply(expr.right(), expr.left().derived());
        }
    };



    template<class V, class Op2, typename T, class Traits, int Backend>
    class Eval<Binary<
                    Tensor<V, 1>,
                    Binary<Number<T>, Unary<Tensor<V, 1>, Op2>, Multiplies>,
                    Minus
                    >, Traits, Backend> {
    public:
        using Expr = Binary<
                            Tensor<V, 1>,
                            Binary<Number<T>, Unary<Tensor<V, 1>, Op2>, Multiplies>,
                            Minus
                        >;

        using Scalar = typename Traits::Scalar;
        using Result = V;

        static Result apply(const Expr &expr)
        {
            Result result;
            apply(expr, result);
            return result;
        }

        static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            // result = (v1) op1 (alpha * unary(v2, op2));
            auto &&v1 = Eval<Tensor<V, 1>, Traits>::apply(expr.left());
            auto &&v2 = Eval<Tensor<V, 1>, Traits>::apply(expr.right().right().expr());

            const Minus &op1 = expr.operation();
            const Op2 &op2 = expr.right().right().operation();
            const Scalar alpha = expr.right().left();

            if(result.is_alias(v1)) {
                assert(false && "avoid using the same object for result and v1");
                V temp;
                temp = v2;
                temp.transform(op2);
                temp.scale(alpha);
                apply_aux(result, temp, op1, result);

            } if(result.is_alias(v2)) {
                result.transform(op2);
                result.scale(alpha);
                apply_aux(v1, result, op1, result);
            } else {
                result = v2;
                result.transform(op2);
                result.scale(alpha);
                apply_aux(v1, result, op1, result);
            }

            UTOPIA_TRACE_END(expr);
        }

        template<class MLeft, class MRight, class MResult>
        inline static void apply_aux(MLeft &&left, MRight &&right, const Minus &op, MResult &result)
        {
            EvalBinaryAux<MResult>::apply(std::forward<MLeft>(left), std::forward<MRight>(right), op, result);
        }

    };


    template<class V, class Op2, typename T, class Traits, int Backend>
    class Eval< Assign<
                    Tensor<V, 1>,
                    Binary<
                        Tensor<V, 1>,
                        Binary<Number<T>, Unary<Tensor<V, 1>, Op2>, Multiplies>,
                        Minus
                    >
                    >, Traits, Backend> {
    public:
        using RightExpr = utopia::Binary<
                            Tensor<V, 1>,
                            Binary<Number<T>, Unary<Tensor<V, 1>, Op2>, Multiplies>,
                            Minus>;

        using Expr = utopia::Assign<Tensor<V, 1>, RightExpr>;
        using Scalar = typename Traits::Scalar;

        static void apply(const Expr &expr) {
            Eval<RightExpr, Traits, Backend>::apply(expr.right(), expr.left().derived());
        }

    };


    //FIXME find better way to avoid specifying all the cases

    template<class T, class Traits, int Backend>
    class Eval<InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Plus>, Plus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Plus>, Plus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&result = Eval<Tensor<T, 1>, Traits>::apply(expr.left());
            result = expr.left() + expr.right().left() + expr.right().right();
            UTOPIA_TRACE_END(expr);
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Minus>, Minus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Minus>, Minus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&result = Eval<Tensor<T, 1>, Traits>::apply(expr.left());
            result = expr.left() - expr.right().left() - expr.right().right();
            UTOPIA_TRACE_END(expr);
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Minus>, Plus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Minus>, Plus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&result = Eval<Tensor<T, 1>, Traits>::apply(expr.left());
            result = expr.left() + expr.right().left() - expr.right().right();
            UTOPIA_TRACE_END(expr);
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Plus>, Minus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<T, 1>, Binary<Tensor<T, 1>, Tensor<T, 1>, Plus>, Minus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&result = Eval<Tensor<T, 1>, Traits>::apply(expr.left());
            result = expr.left() - expr.right().left() + expr.right().right();
            UTOPIA_TRACE_END(expr);
        }
    };


}

#endif //UTOPIA_UTOPIA_EVAL_TERNARY_HPP
