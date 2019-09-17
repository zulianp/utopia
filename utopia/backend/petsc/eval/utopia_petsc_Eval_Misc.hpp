#ifndef UTOPIA_PETSC_EVAL_MISC_HPP
#define UTOPIA_PETSC_EVAL_MISC_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {

        template<class Left, class Right, typename ScalarT, class Traits>
        class Eval<Binary<
                            Binary<Number<ScalarT>, Left,  Multiplies>,
                            Binary<Number<ScalarT>, Right, Multiplies>,
                            Plus
                         >,
                         Traits,
                         PETSC> {
        public:
            typedef utopia::Binary<
                                Binary<Number<ScalarT>, Left,  Multiplies>,
                                Binary<Number<ScalarT>, Right, Multiplies>,
                                Plus> Expr;

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
            {
                EXPR_TYPE(Traits, Expr) result = Eval<Right, Traits>::apply(expr.right().right());

                UTOPIA_TRACE_BEGIN(expr);

                result.axpby(
                    expr.left().left(),
                    Eval<Left, Traits>::apply(expr.left().right()),
                    expr.right().left()
                );

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
                Traits,
                PETSC> {
        public:
            typedef utopia::Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus> Expr;

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
            {
                EXPR_TYPE(Traits, Expr) result;

                UTOPIA_TRACE_BEGIN(expr);

                Eval<Tensor<M, 2>, Traits>::apply(expr.left().left()).multiply_add(
                    Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right()),
                    Eval<Tensor<V2, 1>, Traits>::apply(expr.right()),
                    result
                );

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus>,
                Traits,
                PETSC> {
        public:
            typedef utopia::Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus> Expr;

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
            {
                EXPR_TYPE(Traits, Expr) result;

                UTOPIA_TRACE_BEGIN(expr);

                // UTOPIA_BACKEND(Traits).mat_mult_add(
                //     result,
                //     Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left()),
                //     Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
                //     Eval<Tensor<V1, 1>, Traits>::apply(expr.left())
                // );


                Eval<Tensor<M, 2>, Traits>::apply(expr.right().left()).multiply_add(
                    Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
                    Eval<Tensor<V1, 1>, Traits>::apply(expr.left()),
                    result
                );

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
                Traits,
                PETSC> {
        public:
            typedef utopia::Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus> Expr;

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
            {
                EXPR_TYPE(Traits, Expr) result;

                UTOPIA_TRACE_BEGIN(expr);

                // UTOPIA_BACKEND(Traits).mat_mult_t_add(
                //     result,
                //     Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left().expr()),
                //     Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right()),
                //     Eval<Tensor<V2, 1>, Traits>::apply(expr.right())
                // );

                Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left().expr()).multiply_add(
                    Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right()),
                    Eval<Tensor<V2, 1>, Traits>::apply(expr.right()),
                    result
                );

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus>,
                Traits,
                PETSC> {
        public:
            typedef utopia::Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus> Expr;

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
            {
                EXPR_TYPE(Traits, Expr) result;

                UTOPIA_TRACE_BEGIN(expr);

                // UTOPIA_BACKEND(Traits).mat_mult_t_add(
                //     result,
                //     Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left().expr()),
                //     Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
                //     Eval<Tensor<V1, 1>, Traits>::apply(expr.left())
                // );

                Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left().expr()).multiply_add(
                    Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
                    Eval<Tensor<V1, 1>, Traits>::apply(expr.left()),
                    result
                );

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };
}

#endif //UTOPIA_PETSC_EVAL_MISC_HPP
