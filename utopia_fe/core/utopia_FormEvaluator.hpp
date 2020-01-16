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
        template<class Form, class Derived, int Order>
        static void eval_bilinear(const Form &expr, Tensor<Derived, Order> &tensor, const bool reset_tensor)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_bilinear(expr);
            ctx.init_tensor(tensor.derived(), reset_tensor);
            FormEval<Form, BAKEND_FLAG>::apply(expr, tensor.derived(), ctx);
        }

        template<class Form, class Derived, int Order>
        static void eval_linear(const Form &expr, Tensor<Derived, Order> &tensor, const bool reset_tensor)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_linear(expr);
            eval(expr, tensor.derived(), ctx, reset_tensor);
        }

        template<class Left, class Right, class Matrix, class Vector>
        static void eval_equation(const Equality<Left, Right> &expr, Tensor<Matrix, 2> &mat, Tensor<Vector, 1> &vec, const bool reset_tensors)
        {
            AssemblyContext<BAKEND_FLAG> ctx;
            ctx.init_bilinear(expr.left());
            ctx.init_tensor(mat.derived(), reset_tensors);
            ctx.init_tensor(vec.derived(), reset_tensors);
            FormEval<Equality<Left, Right>, BAKEND_FLAG>::apply(expr, mat.derived(), vec.derived(), ctx);
        }

        template<class Form, class Derived, int Order>
        static void eval(const Form &expr, Tensor<Derived, Order> &tensor, AssemblyContext<BAKEND_FLAG> &ctx, const bool reset_tensor)
        {
            ctx.init_tensor(tensor.derived(), reset_tensor);
            FormEval<Form, BAKEND_FLAG>::apply(expr, tensor.derived(), ctx);
        }

        template<class Form, class Derived, int Order>
        static void eval(Form &expr, Tensor<Derived, Order> &tensor, AssemblyContext<BAKEND_FLAG> &ctx, const bool reset_tensor)
        {
            ctx.init_tensor(tensor.derived(), reset_tensor);
            FormEval<Form, BAKEND_FLAG>::apply(expr, tensor.derived(), ctx);
        }

        template<class Expr, class Matrix, class Vector>
        static void eval(
            const Expr &expr,
            Tensor<Matrix, 2> &mat,
            Tensor<Vector, 1> &vec,
            AssemblyContext<BAKEND_FLAG> &ctx)
        {
            FormEval<Expr, BAKEND_FLAG>::apply(expr, mat.derived(), vec.derived(), ctx);
        }

        template<class Expr, class Matrix, class Vector>
        static void eval(
            Expr &expr,
            Tensor<Matrix, 2> &mat,
            Tensor<Vector, 1> &vec,
            AssemblyContext<BAKEND_FLAG> &ctx)
        {
            FormEval<Expr, BAKEND_FLAG>::apply(expr, mat.derived(), vec.derived(), ctx);
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
