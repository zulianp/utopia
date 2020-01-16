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
            using Expr = utopia::Binary<
                                Binary<Number<ScalarT>, Left,  Multiplies>,
                                Binary<Number<ScalarT>, Right, Multiplies>,
                                Plus>;

            using Result = EXPR_TYPE(Traits, Expr);
            using Scalar = typename Traits::Scalar;

            UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)
            
            // result || right = alpha * left + beta * right
            inline static void apply(const Expr &expr, Result &result)
            {
                UTOPIA_TRACE_BEGIN(expr);

                const Scalar alpha = expr.left().left();
                const Scalar beta  = expr.right().left();

                auto &&right = Eval<Right, Traits>::apply(expr.right().right());
                auto &&left  = Eval<Left, Traits>::apply(expr.left().right());

                if(result.is_alias(right)) {
                    result.axpby(
                        alpha,
                        left,
                        beta
                    );
                } else if(result.is_alias(left)) {
                    result.axpby(
                        beta,
                        right,
                        alpha
                    );

                } else {
                    result.construct(right);
                    result.axpby(
                        alpha,
                        left,
                        beta
                    );
                }

                UTOPIA_TRACE_END(expr);
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
                Traits,
                PETSC> {
        public:
            using Expr = utopia::Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>;
            using Result = EXPR_TYPE(Traits, Expr);

            UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

            inline static void apply(const Expr &expr, Result &result)
            {
                UTOPIA_TRACE_BEGIN(expr);

                auto &&v1 = Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right());
                auto &&v2 = Eval<Tensor<V2, 1>, Traits>::apply(expr.right());

                if(result.is_alias(v1) || result.is_alias(v2)) {

                    Result temp;

                    Eval<Tensor<M, 2>, Traits>::apply(expr.left().left()).multiply_add(
                        v1,
                        v2,
                        temp
                    );

                    result.construct(std::move(temp));

                } else {
                    Eval<Tensor<M, 2>, Traits>::apply(expr.left().left()).multiply_add(
                        v1,
                        v2,
                        result
                    );
                }

                UTOPIA_TRACE_END(expr);
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus>,
                Traits,
                PETSC> {
        public:
            using Expr = utopia::Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus>;
            using Result = EXPR_TYPE(Traits, Expr);

            UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

            inline static void apply(const Expr &expr, Result &result)
            {
                UTOPIA_TRACE_BEGIN(expr);
                auto &&v2 = Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right());
                auto &&v1 = Eval<Tensor<V1, 1>, Traits>::apply(expr.left());

                if(result.is_alias(v1) || result.is_alias(v2)) {
                    Result temp;

                    Eval<Tensor<M, 2>, Traits>::apply(expr.right().left()).multiply_add(
                        v2,
                        v1,
                        result
                    );

                    result.construct(std::move(temp));

                } else {
                    Eval<Tensor<M, 2>, Traits>::apply(expr.right().left()).multiply_add(
                        v2,
                        v1,
                        result
                    );
                }

                UTOPIA_TRACE_END(expr);
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
                Traits,
                PETSC> {
        public:
            using Expr = utopia::Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>;
            using Result = EXPR_TYPE(Traits, Expr);

            UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

            inline static void apply(const Expr &expr, Result &result)
            {
                UTOPIA_TRACE_BEGIN(expr);

                auto &&v1 = Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right());
                auto &&v2 = Eval<Tensor<V2, 1>, Traits>::apply(expr.right());

                if(result.is_alias(v1) || result.is_alias(v2)) {
                    Result temp;
                    
                    Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left().expr()).multiply_add(
                        v1,
                        v2,
                        temp
                    );

                    result.assign(std::move(temp));

                } else {
                    Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left().expr()).multiply_add(
                        v1,
                        v2,
                        result
                    );
                }

                UTOPIA_TRACE_END(expr);
            }
        };

        template<class M, class V1, class V2, class Traits>
        class Eval<
                Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus>,
                Traits,
                PETSC> {
        public:
            using Expr = utopia::Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus>;
            using Result = EXPR_TYPE(Traits, Expr);

            UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr, Result &result)
            {
                UTOPIA_TRACE_BEGIN(expr);

                auto &&v2 = Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right());
                auto &&v1 = Eval<Tensor<V1, 1>, Traits>::apply(expr.left());

                if(result.is_alias(v1) || result.is_alias(v2)) {

                    Result temp;
                    Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left().expr()).multiply_add(
                        v2,
                        v1,
                        temp
                    );

                    result.assign(std::move(temp));

                } else {
                    Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left().expr()).multiply_add(
                        v2,
                        v1,
                        result
                    );
                }

                UTOPIA_TRACE_END(expr);
                return result;
            }
        };
}

#endif //UTOPIA_PETSC_EVAL_MISC_HPP
