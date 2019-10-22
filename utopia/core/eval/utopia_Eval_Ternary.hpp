#ifndef UTOPIA_UTOPIA_EVAL_TERNARY_HPP
#define UTOPIA_UTOPIA_EVAL_TERNARY_HPP

#include "utopia_Eval_Binary.hpp"


namespace utopia {

    template<class T1, class T2, class T3, int Order, class Op1, class Op2>
    using TernaryExpr = utopia::Binary<utopia::Binary<Tensor<T1, Order>, Tensor<T2, Order>, Op1>, Tensor<T3, Order>, Op2>;

    template<class Result, class T1, class T2, class T3, int Order, class Op1, class Op2, class Traits, int Backend>
    class Eval<Assign<Tensor<Result, Order>, TernaryExpr<T1, T2, T3, Order, Op1, Op2>>, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Tensor<Result, Order>, TernaryExpr<T1, T2, T3, Order, Op1, Op2>>;
        
        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &result = Eval<Tensor<Result, Order>, Traits>::apply(expr.left());
            auto &&t1    = Eval<Tensor<T1, Order>, Traits>::apply(expr.right().left().left());
            auto &&t2    = Eval<Tensor<T1, Order>, Traits>::apply(expr.right().left().right());
            auto &&t3    = Eval<Tensor<T3, Order>, Traits>::apply(expr.right().right());

            const auto op1 = expr.right().left().operation();
            const auto op2 = expr.right().operation();

            const static bool order_matters = !eval_order_changable<Op1, Op2>::value;

            if(result.same_object(t3)) {
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

            UTOPIA_TRACE_END_SPECIALIZED(expr);
        }

    private:
        template<class AuxR, class AuxT1, class AuxT2, class AuxT3, class AuxOp1, class AuxOp2>
        inline static void apply_aux(AuxR &result, AuxT1 &&t1, AuxT2 &&t2, AuxT3 &&t3, const AuxOp1 &op1, const AuxOp2 &op2)
        {
            EvalBinaryAux<AuxR>::apply(t1, t2, op1, result);
            EvalBinaryAux<AuxR>::apply(result, t3, op2, result);
        }
    };


    //result = a op1 (alpha * b) op2 (beta * c);
    //Assign<Vec, Minus<Plus<Vec, Multiplies<Number, Vec>>, Multiplies<Number, Vec>>>
    template<class V, typename T, class Op1, class Op2>
    using ScaledTernaryExpr = utopia::Assign<
                                Tensor<V, 1>, 
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
                                    Op2>
                                >;

    template<class V, typename T, class Traits, int Backend>
    class Eval< ScaledTernaryExpr<V, T, Plus, Minus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<V, 1>;
        using Scalar = typename Traits::Scalar;

        static void apply(const ScaledTernaryExpr<V, T, Plus, Minus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            //result = (a + (alpha * b)) - (beta * c);
            auto &result = Eval<T1, Traits>::apply(expr.left());
            auto &a = Eval<T1, Traits>::apply(expr.right().left().left());
            const Scalar alpha = expr.right().left().right().left();
            auto &b = Eval<T1, Traits>::apply(expr.right().left().right().right());
            const Scalar beta = expr.right().right().left(); 
            auto &c = Eval<T1, Traits>::apply(expr.right().right().right());

            if(result.same_object(b)) {
                result.scale(alpha);
                result.axpy(1.0, a);
                result.axpy(-beta, c);
            } else if(result.same_object(c)) {
                result.scale(-beta);
                result.axpy(alpha, b);
                result.axpy(1.0, a);
            } else if(result.same_object(a)) {
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
    class Eval< ScaledTernaryExpr<V, T, Minus, Plus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<V, 1>;
        using Scalar = typename Traits::Scalar;

        static void apply(const ScaledTernaryExpr<V, T, Minus, Plus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            //result = (a - (alpha * b)) + (beta * c);
            auto &result = Eval<T1, Traits>::apply(expr.left());
            auto &a = Eval<T1, Traits>::apply(expr.right().left().left());
            const Scalar alpha = expr.right().left().right().left();
            auto &b = Eval<T1, Traits>::apply(expr.right().left().right().right());
            const Scalar beta = expr.right().right().left(); 
            auto &c = Eval<T1, Traits>::apply(expr.right().right().right());

            if(result.same_object(b)) {
                result.scale(-alpha);
                result.axpy(1.0, a);
                result.axpy(beta, c);
            } else if(result.same_object(c)) {
                result.scale(beta);
                result.axpy(-alpha, b);
                result.axpy(1.0, a);
            } else if(result.same_object(a)) {
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


    
    //Assign<PetscMatrix, Plus<PetscMatrix, Multiply<PetscMatrix, PetscMatrix>>>
    template<class M, class Op, class Traits, int Backend>
    class Eval< Assign<M, Binary<M, Multiply<M, M>, Op>>, Traits, Backend> {
    public:
        static void apply(const Assign<M, Binary<M, Multiply<M, M>, Op>> &expr)
        {
            // result = A op (B * C)
            auto &result = Eval<M, Traits>::apply(expr.left());
            auto &&A = Eval<M, Traits>::apply(expr.right().left());
            auto &&B = Eval<M, Traits>::apply(expr.right().right().left());
            auto &&C = Eval<M, Traits>::apply(expr.right().right().right());
            auto && op = expr.right().operation();

            if(result.same_object(A) || result.same_object(B) || result.same_object(C)) {
                typename Unwrap<M>::Type temp;
                B.multiply(C, temp);
                apply_aux(A, std::move(temp), op, result);
            } else {
                B.multiply(C, result);
                apply_aux(A, result, op, result);
            }
        }

        template<class MLeft, class MRight, class MResult>
        inline static void apply_aux(const MLeft &left, MRight &&right, const Op &op, MResult &result)
        {
            EvalBinaryAux<MResult>::apply(left, std::forward<MRight>(right), op, result);
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
        using Expr = Assign<
                        Tensor<V, 1>,
                        Binary<
                            Tensor<V, 1>,
                            Binary<Number<T>, Unary<Tensor<V, 1>, Op2>, Multiplies>,
                            Minus
                        >
                    >;
        using Scalar = typename Traits::Scalar;

        static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            // result = (v1) op1 (alpha * unary(v2, op2)); 
            auto &result = Eval<Tensor<V, 1>, Traits>::apply(expr.left());
            auto &&v1 = Eval<Tensor<V, 1>, Traits>::apply(expr.right().left());
            auto &&v2 = Eval<Tensor<V, 1>, Traits>::apply(expr.right().right().right().expr());

            const Minus &op1 = expr.right().operation();
            const Op2 &op2 = expr.right().right().right().operation();
            const Scalar alpha = expr.right().right().left();

            if(result.same_object(v1)) {
                assert(false && "avoid using the same object for result and v1");
                V temp;
                temp = v2;
                temp.transform(op2);
                temp.scale(alpha);
                apply_aux(result, temp, op1, result);

            } if(result.same_object(v2)) {
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

}

#endif //UTOPIA_UTOPIA_EVAL_TERNARY_HPP
