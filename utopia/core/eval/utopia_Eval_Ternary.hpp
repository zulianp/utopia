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

}

#endif //UTOPIA_UTOPIA_EVAL_TERNARY_HPP
