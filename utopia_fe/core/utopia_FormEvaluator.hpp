#ifndef UTOPIA_FE_FORM_EVALUATOR_HPP
#define UTOPIA_FE_FORM_EVALUATOR_HPP

#include "utopia_FormExpressions.hpp"
#include "utopia_Traits.hpp"
#include "utopia_fe_lang.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FormEval.hpp"

namespace utopia {

    // Compile-time member detection flags.
    class FeatureYes { public: char c[1]; };
    class FeatureNo  { public: char c[2]; };

    template<class T>
    FeatureYes DetectIsFE(decltype(&T::is_fe));

    template<class T>
    FeatureNo DetectIsFE(...);

    //, sizeof(DetectIsFE<Form>(0))

    template<int BAKEND_FLAG>
    class FormEvaluator {
    public:
        template<class Form, class Tensor, int Order>
        static void eval_bilinear(const Form &expr, Wrapper<Tensor, Order> &tensor, const bool reset_tensor)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_bilinear(expr);
            ctx.init_tensor(tensor, reset_tensor);
            FormEval<Form, BAKEND_FLAG>::apply(expr, tensor, ctx);
        }

        template<class Form, class Tensor, int Order>
        static void eval_linear(const Form &expr, Wrapper<Tensor, Order> &tensor, const bool reset_tensor)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_linear(expr);
            eval(expr, tensor, ctx, reset_tensor);
        }

        template<class Left, class Right, class Matrix, class Vector>
        static void eval_equation(const Equality<Left, Right> &expr, Wrapper<Matrix, 2> &mat, Wrapper<Vector, 1> &vec, const bool reset_tensors)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_bilinear(expr.left());
            ctx.init_tensor(mat, reset_tensors);
            ctx.init_tensor(vec, reset_tensors);
            FormEval<Equality<Left, Right>, BAKEND_FLAG>::apply(expr, mat, vec, ctx);
        }

        template<class Form, class Tensor, int Order>
        static void eval(const Form &expr, Wrapper<Tensor, Order> &tensor, AssemblyContext<BAKEND_FLAG> &ctx, const bool reset_tensor)
        {
            ctx.init_tensor(tensor, reset_tensor);
            FormEval<Form, BAKEND_FLAG>::apply(expr, tensor, ctx);
        }


        template<class Expr, class Matrix, class Vector>
        static void eval(
            const Expr &expr,
            Wrapper<Matrix, 2> &mat,
            Wrapper<Vector, 1> &vec,
            AssemblyContext<BAKEND_FLAG> &ctx)
        {
            FormEval<Expr, BAKEND_FLAG>::apply(expr, mat, vec, ctx);
        }



        template<class Expr, typename T>
        static void eval(
            const Expr &expr,
            Number<T> &value,
            AssemblyContext<BAKEND_FLAG> &ctx)
        {
            FormEval<Expr, BAKEND_FLAG>::apply(expr, value, ctx);
        }
    };
}

#endif //UTOPIA_FE_FORM_EVALUATOR_HPP
