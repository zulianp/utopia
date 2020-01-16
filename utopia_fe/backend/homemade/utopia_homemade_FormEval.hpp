// #ifndef UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP
// #define UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP

// #include "utopia_Base.hpp"
// #include "utopia_homemade_FEForwardDeclarations.hpp"

// #include "utopia_FormEval.hpp"
// #include "utopia_FEEval.hpp"

// #include "utopia_homemade_AssemblyContext.hpp"
// #include "utopia_homemade_FunctionSpace.hpp"


// namespace utopia {
//     template<class Form>
//     class FormEval<Form, HOMEMADE> {
//     public:
//         typedef utopia::Traits<HMFESpace> Traits;
//         FormEval() { }

//         template<class Expr, class Tensor>
//         static void apply(
//                     const Integral<Expr> &expr,
//                     Tensor &t,
//                     AssemblyContext<HOMEMADE> &ctx)
//         {
//             if(expr.has_block_id() && ctx.block_id() != expr.block_id()) {
//                 return;
//             }

//             auto &&r = FEEval<Integral<Expr>, Traits, HOMEMADE, QUAD_DATA_NO>::apply(expr, ctx);
//             t = r;
//         }


//         template<class Left, class Right, class Matrix, class Vector>
//         static void apply(
//                     const Equality<Left, Right> &expr,
//                     Wrapper<Matrix, 2> &mat,
//                     Wrapper<Vector, 1> &vec,
//                     AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.left(),  mat, ctx);
//             apply(expr.right(), vec, ctx);
//         }


//         template<class Left, class Right, class Tensor>
//         static void apply(
//             const Binary<Left, Right, Plus> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.left(),  result, ctx);

//             Tensor right = zeros(size(result));
//             apply(expr.right(), result, ctx);
//             result += right;
//         }

//         template<class Left, class Right, class Tensor>
//         static void apply(
//             const Binary<Number<Left>, Right, Multiplies> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.right(), result, ctx);
//             result *= expr.left();
//         }

//         template<class Left, class Right, class Tensor>
//         static void apply(
//             const Binary<Left, Right, Minus> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.left(), result, ctx);

//             Tensor right = zeros(size(result));
//             apply(expr.right(), result, ctx);
//             result -= right;
//         }

//         template<class Left, class Right, class Op, class Tensor>
//         static void apply(
//             const Binary<Left, Right, Op> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.left(), result, ctx);

//             Tensor right = zeros(size(result));
//             apply(expr.right(), result, ctx);

//             result = expr.operation().apply(result, right);
//         }

//         template<class Expr, class Tensor>
//         static void apply(
//             const Unary<Expr, Minus> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.expr(), result, ctx);
//             result = -result;
//         }

//         template<class Expr, class Tensor>
//         static void apply(
//             const Negate<Expr> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.expr(), result, ctx);
//             result = -result;
//         }

//         template<class Expr, class Tensor>
//         static void apply(
//             const Unary<Expr, Abs> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.expr(), result, ctx);
//             result = abs(result);
//         }

//         template<class Expr, class Tensor>
//         static void apply(
//             const Unary<Expr, Sqrt> &expr,
//             Tensor &result,
//             AssemblyContext<HOMEMADE> &ctx)
//         {
//             apply(expr.expr(), result, ctx);
//             result = sqrt(result);
//         }
//     };
// }

// #endif //UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP
