//TODO
#ifndef UTOPIA_EVAL_PETSC_HPP
#define UTOPIA_EVAL_PETSC_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Eval_Rename.hpp"
#include "utopia_petsc_Eval_KroneckerProduct.hpp"
// #include "utopia_petsc_Eval_NZZXRow.hpp"
// #include "utopia_petsc_Eval_Inverse.hpp"
// #include "utopia_petsc_Eval_Factory.hpp"
// #include "utopia_petsc_Eval_Residual.hpp"
// #include "utopia_petsc_Eval_DotOpDot.hpp"
// #include "utopia_petsc_Eval_Blocks.hpp"
// #include "utopia_petsc_EvalDotVecVecs.hpp"
// #include "utopia_petsc_EvalMatGetCol.hpp"
// #include "utopia_petsc_Eval_VecUniqueSortSerial.hpp"
// #include "utopia_petsc_Eval_Chop.hpp"

// #ifdef WITH_SLEPC
// #include "utopia_petsc_Eval_Cond.hpp"
// #endif //WITH_SLEPC

// /*! @file
// * Petsc language extensions
// */

// namespace utopia {
//     template<class Left, class Right, class Traits>
//     class Eval<Construct<Left, LocalDiagBlock<Right> >, Traits, PETSC> {
//     public:

//         inline static void apply(const Construct<Left, LocalDiagBlock<Right> > & expr)
//         {
//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).build_local_diag_block(
//                 Eval<Left,  Traits>::apply(expr.left()),
//                 Eval<Right, Traits>::apply(expr.right().expr())
//                 );

//             UTOPIA_TRACE_END(expr);
//         }
//     };

//     // mat-mat-mat multiplication
//     template<class M1, class M2, class M3, class Traits>
//     class Eval< Multiply< Multiply< Wrapper<M1, 2>, Wrapper<M2, 2> >, Wrapper<M3, 2> >, Traits, PETSC> {
//         public:
//             typedef utopia::Multiply< Multiply< Wrapper<M1, 2>, Wrapper<M2, 2> >, Wrapper<M3, 2> > Expr;

//             inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//             {
//                 EXPR_TYPE(Traits, Expr) result;

//                 UTOPIA_TRACE_BEGIN(expr);

//                 //Performs optimal triple product
//                 UTOPIA_BACKEND(Traits).triple_product(
//                     result,
//                     Eval<Wrapper<M1, 2>, Traits>::apply(expr.left().left()),
//                     Eval<Wrapper<M2, 2>, Traits>::apply(expr.left().right()),
//                     Eval<Wrapper<M3, 2>, Traits>::apply(expr.right())
//                     );

//                 UTOPIA_TRACE_END(expr);
//                 // assert(result.same_type(Eval<Wrapper<M3, 2>, Traits>::apply(expr.right())));
//                 return result;
//             }
//     };

//     //! [pattern matching and optimizations]

//     /*!
//     * @brief Triple product (m1^T * m2 * m1 := transpose(m1) * m2 * m1) optimization for the petsc backend
//     */
//     template<class M1, class M2, class Traits>
//     class Eval<
//         //The pattern to match
//         Multiply< Multiply<Transposed<M1>, M2>, M1>,
//         //Type information
//         Traits,
//         //Restriction to the backend with the PETSC tag.
//         PETSC> {
//     public:
//         typedef utopia::Multiply< Multiply<Transposed<M1>, M2>, M1> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             //Check if left and right operands are the same object
//             if(&expr.left().left().expr() == &expr.right()) {
//                 //Perform optimal triple product
//                 UTOPIA_BACKEND(Traits).triple_product_ptap(
//                     result,
//                     Eval<M1, Traits>::apply(expr.left().right()),
//                     Eval<M2, Traits>::apply(expr.right())
//                     );

//             } else {
//                 //Perform general triple product
//                 UTOPIA_BACKEND(Traits).triple_product(
//                     result,
//                     Eval<Transposed<M1>, Traits>::apply(expr.left().left()),
//                     Eval<M2, Traits>::apply(expr.left().right()),
//                     Eval<M1, Traits>::apply(expr.right())
//                     );
//             }

//             // assert( result.same_type(Eval<M1, Traits>::apply(expr.right())) );
//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };



//     /*!
//     * @brief Triple product (m1 * m2 * m1^T := m1 * m2 * transpose(m1)) optimization for the petsc backend
//     */
//     template<class M1, class M2, class Traits>
//     class Eval<
//         //The pattern to match
//         Multiply< Multiply<M1, M2>, Transposed<M1>>,
//         //Type information
//         Traits,
//         //Restriction to the backend with the PETSC tag.
//         PETSC> {
//     public:
//         typedef utopia::Multiply< Multiply<M1, M2>, Transposed<M1>> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             //Check if left and right operands are the same object
//             if(&expr.left().left() == &expr.right().expr()) {
//                 //Perform optimal triple product
//                 UTOPIA_BACKEND(Traits).triple_product_rart(
//                     result,
//                     Eval<M1, Traits>::apply(expr.left().right()),
//                     Eval<M2, Traits>::apply(expr.right().expr())
//                 );

//             } else {
//                 //Perform general triple product
//                 UTOPIA_BACKEND(Traits).triple_product(
//                     result,
//                     Eval<M1, Traits>::apply(expr.left().left()),
//                     Eval<M2, Traits>::apply(expr.left().right()),
//                     Eval<Transposed<M1>, Traits>::apply(expr.right())
//                 );
//             }

//             // assert( result.same_type(Eval<M1, Traits>::apply(expr.right())) );
//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

//     //! [pattern matching and optimizations]

//     template<class Left, class Right, class Traits>
//     class Eval<LocalRedistribute<Left, Right>, Traits, PETSC> {
//     public:

//         inline static EXPR_TYPE(Traits, Left) apply(const LocalRedistribute<Left, Right> &expr)
//         {
//             EXPR_TYPE(Traits, Left) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).build_local_redistribute(
//                 result,
//                 Eval<Left,  Traits>::apply(expr.left()),
//                 Eval<Right, Traits>::apply(expr.right())
//                 );


//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };



//     template<class Left, class Right, typename ScalarT, class Traits>
//     class Eval<Binary<
//                         Binary<Number<ScalarT>, Left,  Multiplies>,
//                         Binary<Number<ScalarT>, Right, Multiplies>,
//                         Plus
//                      >,
//                      Traits,
//                      PETSC> {
//     public:
//         typedef utopia::Binary<
//                             Binary<Number<ScalarT>, Left,  Multiplies>,
//                             Binary<Number<ScalarT>, Right, Multiplies>,
//                             Plus> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result = Eval<Right, Traits>::apply(expr.right().right());

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).axpby(
//                     result,
//                     expr.left().left(),
//                     Eval<Left,  Traits>::apply(expr.left().right()),
//                     expr.right().left()
//             );

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

//     template<class M, class V1, class V2, class Traits>
//     class Eval<
//             Binary<Multiply<Wrapper<M, 2>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus>,
//             Traits,
//             PETSC> {
//     public:
//         typedef utopia::Binary<Multiply<Wrapper<M, 2>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).mat_mult_add(
//                 result,
//                 Eval<Wrapper<M, 2>,  Traits>::apply(expr.left().left()),
//                 Eval<Wrapper<V1, 1>, Traits>::apply(expr.left().right()),
//                 Eval<Wrapper<V2, 1>, Traits>::apply(expr.right())
//             );

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

//     template<class M, class V1, class V2, class Traits>
//     class Eval<
//             Binary<Wrapper<V1, 1>, Multiply<Wrapper<M, 2>, Wrapper<V2, 1>>, Plus>,
//             Traits,
//             PETSC> {
//     public:
//         typedef utopia::Binary<Wrapper<V1, 1>, Multiply<Wrapper<M, 2>, Wrapper<V2, 1>>, Plus> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).mat_mult_add(
//                 result,
//                 Eval<Wrapper<M, 2>,  Traits>::apply(expr.right().left()),
//                 Eval<Wrapper<V2, 1>, Traits>::apply(expr.right().right()),
//                 Eval<Wrapper<V1, 1>, Traits>::apply(expr.left())
//             );

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

//     template<class M, class V1, class V2, class Traits>
//     class Eval<
//             Binary<Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus>,
//             Traits,
//             PETSC> {
//     public:
//         typedef utopia::Binary<Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).mat_mult_t_add(
//                 result,
//                 Eval<Wrapper<M, 2>,  Traits>::apply(expr.left().left().expr()),
//                 Eval<Wrapper<V1, 1>, Traits>::apply(expr.left().right()),
//                 Eval<Wrapper<V2, 1>, Traits>::apply(expr.right())
//             );

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };

//     template<class M, class V1, class V2, class Traits>
//     class Eval<
//             Binary<Wrapper<V1, 1>, Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V2, 1>>, Plus>,
//             Traits,
//             PETSC> {
//     public:
//         typedef utopia::Binary<Wrapper<V1, 1>, Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V2, 1>>, Plus> Expr;

//         inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
//         {
//             EXPR_TYPE(Traits, Expr) result;

//             UTOPIA_TRACE_BEGIN(expr);

//             UTOPIA_BACKEND(Traits).mat_mult_t_add(
//                 result,
//                 Eval<Wrapper<M, 2>,  Traits>::apply(expr.right().left().expr()),
//                 Eval<Wrapper<V2, 1>, Traits>::apply(expr.right().right()),
//                 Eval<Wrapper<V1, 1>, Traits>::apply(expr.left())
//             );

//             UTOPIA_TRACE_END(expr);
//             return result;
//         }
//     };
// }

#endif //UTOPIA_EVAL_PETSC_HPP
