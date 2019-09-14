//TODO
// #ifndef UTOPIA_BLAS_EVAL_MULTIPLY_HPP
// #define UTOPIA_BLAS_EVAL_MULTIPLY_HPP

// #include "utopia_Eval_Empty.hpp"

// namespace utopia {

//     //C := alpha * op( A ) * op( B )
//     template<typename Alpha, class A, class B, class Traits>
//     class Eval<
//             Multiply<
//                     Binary<Number<Alpha>, A, Multiplies>,
//                     Tensor<B, 2>
//                 >,
//             Traits, utopia::BLAS
//             > {
//     public:
//         typedef utopia::Multiply< Binary<Number<Alpha>, A, Multiplies>, Tensor<B, 2>> Expr;
//         typedef typename TypeAndFill<Traits, Multiply<A, Tensor<B, 2>> >::Type C;

//         inline static C apply(const Expr &expr) {
//             C c;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).gemm(
//                 //C
//                 c,
//                 //beta
//                 0.,
//                 //alpha
//                 expr.left().left(),
//                 //transpose_A
//                 false,
//                 //A
//                 Eval<A, Traits>::apply(expr.left().right()),
//                 //transpose_B
//                 false,
//                 //B
//                 Eval<Tensor<B, 2>, Traits>::apply(expr.right())
//             );

//             UTOPIA_TRACE_END(expr);
//             return c;
//         }
//     };

//     //C := alpha * op( A ) * op( B ) + beta * C,
//     template<typename Alpha, class A, class B, typename Beta, class C, class Traits>
//     class Eval<
//             Binary<
//                 Multiply<
//                         Binary<Number<Alpha>, A, Multiplies>,
//                         B
//                     >,
//                 Binary<
//                         Number<Beta>,
//                         Tensor<C, 2>,
//                         Multiplies
//                     >,
//                 Plus
//             >,
//             Traits,
//             utopia::BLAS
//             > {
//     public:
//         typedef utopia::Binary<
//                 Multiply<
//                         Binary<Number<Alpha>, A, Multiplies>,
//                         B
//                     >,
//                 Binary<
//                         Number<Beta>,
//                         Tensor<C, 2>,
//                         Multiplies
//                     >,
//                 Plus
//             > Expr;


//         typedef typename TypeAndFill<Traits, Tensor<C, 2>>::Type Result;

//         inline static Result apply(const Expr &expr) {
//             Result result = Eval<Tensor<C, 2>, Traits>::apply(expr.right().right());

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).gemm(
//                 //Result
//                 result,
//                 //beta
//                 expr.right().left(),
//                 //alpha
//                 expr.left().left().left(),
//                 //transpose_A
//                 false,
//                 //A
//                 Eval<A, Traits>::apply(expr.left().left().right()),
//                 //transpose_B
//                 false,
//                 //B
//                 Eval<B, Traits>::apply(expr.left().right())
//             );

//             // std::cout << tree_format(expr.get_class()) << std::endl;

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

// }

// #endif //UTOPIA_BLAS_EVAL_MULTIPLY_HPP
