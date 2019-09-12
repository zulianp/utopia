//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_REDUCE_HPP
#define UTOPIA_UTOPIA_EVAL_REDUCE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Each.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Operations.hpp"

namespace utopia {

    template <class Expr, class Operation, class Traits, int Backend>
    class Eval<Reduce<Expr, Operation>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Reduce<Expr, Operation> &expr)
        {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<Expr, Traits>::apply(expr.expr()).reduce(expr.operation());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Dot<Left, Right>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Dot<Left, Right> &expr)
        {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).dot(
            //         Eval<Left,  Traits>::apply(expr.expr().left()),
            //         Eval<Right, Traits>::apply(expr.expr().right())
            // );

            result = Eval<Left,  Traits>::apply(expr.expr().left()).dot(
                Eval<Right, Traits>::apply(expr.expr().right())
            );


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Norm<Expr, 2>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Norm<Expr, 2> &expr)
        {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<Expr, Traits>::apply(expr.expr()).norm2();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Norm<Expr, 1>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Norm<Expr, 1> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).norm1(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Norm<Expr, INFINITY_NORM_TAG>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Norm<Expr, INFINITY_NORM_TAG> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).norm_infty(
            //         Eval<Expr, Traits>::apply(expr.expr())
            // );

            result = Eval<Expr, Traits>::apply(expr.expr()).norm_infty();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<ReduceApproxEqual<Left, Right>, Traits, Backend> {
    public:
        inline static bool apply(const ReduceApproxEqual<Left, Right> &expr) {
            bool result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).compare(
            //         Eval<Left,  Traits>::apply(expr.expr().left()),
            //         Eval<Right, Traits>::apply(expr.expr().right()),
            //         expr.expr().operation());


            result = 
                Eval<Left,  Traits>::apply(expr.expr().left()).equals(
                Eval<Right, Traits>::apply(expr.expr().right()),
                expr.expr().operation().tol()
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Reduce< Diag<Expr>, Plus>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;
        inline static Scalar apply(const Reduce< Diag<Expr>, Plus> &expr)
        {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).trace(
                    Eval<Expr, Traits>::apply(expr.expr().expr())
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Trace<Expr>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Trace<Expr> &expr)
        {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).trace(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };


    template<class Derived, typename T, class Traits, int Backend>
    class Eval<Reduce<Tensor<Derived, 2>, PlusIsNonZero<T>>, Traits, Backend> {
    public:
        using Expr = utopia::Reduce<Tensor<Derived, 2>, PlusIsNonZero<T>>;
        using Vector = UTOPIA_VECTOR(Derived);

        typedef typename Traits::Scalar Scalar;
        typedef typename Traits::SizeType SizeType;

        inline static SizeType apply(const Expr &expr)
        {
            const auto &op = expr.operation().is_non_zero();
            SizeType result = 0;
            UTOPIA_TRACE_BEGIN(expr);

            each_read(expr.expr().derived(), [&result, &op](const SizeType i, const SizeType j, const Scalar value) {
                UTOPIA_UNUSED(i);
                UTOPIA_UNUSED(j);

                result += op.apply(value);
            });


            Vector v = local_values(1, result);
            result = sum(v);

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_REDUCE_HPP
