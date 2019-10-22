#ifndef UTOPIA_BLAS_EVAL_MULTIPLY_HPP
#define UTOPIA_BLAS_EVAL_MULTIPLY_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    //C := alpha * op( A ) * op( B )
    template<typename Alpha, class A, class B, class Traits>
    class Eval<
            Multiply<
                    Binary<Number<Alpha>, A, Multiplies>,
                    Tensor<B, 2>
                >,
            Traits, utopia::BLAS
            > {
    public:
        typedef utopia::Multiply< Binary<Number<Alpha>, A, Multiplies>, Tensor<B, 2>> Expr;
        typedef typename TypeAndFill<Traits, Tensor<B, 2> >::Type C;

        static_assert(std::is_same<C, B>::value, "output must be a 2nd order tensor");

        inline static C apply(const Expr &expr)
        {
            C c;
            apply(expr, c);
            return c;
        }   

        inline static C apply(const Expr &expr, C &c)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&a = Eval<A, Traits>::apply(expr.left().right());
            auto &&b = Eval<Tensor<B, 2>, Traits>::apply(expr.right());
            const Alpha alpha = expr.left().left();

            // a.gemm(false, alpha, false, b, 0.0, c);
            apply_aux(alpha, a, b, c);

            UTOPIA_TRACE_END(expr);
        }

    private:
        template<class MatT>
        inline static void apply_aux(
            const Number<Alpha> &alpha,
            const Tensor<MatT, 2> &a,
            const Tensor<MatT, 2> &b,
            Tensor<MatT, 2> &c)
        {
            a.derived().gemm(false, alpha, false, b.derived(), 0.0, c.derived());
        }

        template<class MatT>
        inline static void apply_aux(
            const Number<Alpha> &alpha,
            const Number<Alpha> &a,
            const Tensor<MatT, 2> &b,
            Tensor<MatT, 2> &c)
        {
            // a.gemm(false, alpha, false, b, 0.0, c);

            c.derived() = b.derived();
            c.derived() *= (alpha * a);
        }

        //move version
        template<class MatT>
        inline static void apply_aux(
            const Number<Alpha> &alpha,
            const Number<Alpha> &a,
            Tensor<MatT, 2> &&b,
            Tensor<MatT, 2> &c)
        {
            // a.gemm(false, alpha, false, b, 0.0, c);

            c.derived() = std::move(b.derived());
            c.derived() *= (alpha * a);
        }

        template<class MatT>
        inline static void apply_aux(
            const Number<Alpha> &alpha,
            const Tensor<MatT, 2> &a,
            const Number<Alpha> &b,
            Tensor<MatT, 2> &c)
        {
            // a.gemm(false, alpha, false, b, 0.0, c);

            c.derived() = a.derived();
            c.derived() *= (b * alpha);
        }

        //move version
        template<class MatT>
        inline static void apply_aux(
            const Number<Alpha> &alpha,
            Tensor<MatT, 2> &&a,
            const Number<Alpha> &b,
            Tensor<MatT, 2> &c)
        {
            // a.gemm(false, alpha, false, b, 0.0, c);

            c.derived() = std::move(a.derived());
            c.derived() *= (b * alpha);
        }
    };

    //C := alpha * op( A ) * op( B ) + beta * C,
    template<typename Alpha, class A, class B, typename Beta, class C, class Traits>
    class Eval<
            Binary<
                Multiply<
                        Binary<Number<Alpha>, A, Multiplies>,
                        B
                    >,
                Binary<
                        Number<Beta>,
                        Tensor<C, 2>,
                        Multiplies
                    >,
                Plus
            >,
            Traits,
            utopia::BLAS
            > {
    public:
        typedef utopia::Binary<
                Multiply<
                        Binary<Number<Alpha>, A, Multiplies>,
                        B
                    >,
                Binary<
                        Number<Beta>,
                        Tensor<C, 2>,
                        Multiplies
                    >,
                Plus
            > Expr;


        typedef typename TypeAndFill<Traits, Tensor<C, 2>>::Type Result;

        inline static Result apply(const Expr &expr) {
            Result result;
            apply(expr, result);
            return result;
        }

        inline static Result apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&a = Eval<A, Traits>::apply(expr.left().left().right());
            auto &&b = Eval<B, Traits>::apply(expr.left().right());
            auto &&c = Eval<Tensor<C, 2>, Traits>::apply(expr.right().right());

            auto &&alpha = expr.left().left().left();
            auto &&beta  =  expr.right().left();

            if(result.same_object(a)) {
                Result temp = std::move(result);
                result = c;
                temp.gemm(false, alpha, false, b, beta, result);
            } else if(result.same_object(b)) {
                Result temp = std::move(result);
                result = c;
                a.gemm(false, alpha, false, temp, beta, result);
            } else {
                result = c;
                a.gemm(false, alpha, false, b, beta, result);
            }

            UTOPIA_TRACE_END(expr);
        }
    };

}

#endif //UTOPIA_BLAS_EVAL_MULTIPLY_HPP
