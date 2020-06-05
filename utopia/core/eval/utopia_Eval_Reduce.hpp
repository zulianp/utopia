//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_REDUCE_HPP
#define UTOPIA_UTOPIA_EVAL_REDUCE_HPP

//#include "utopia_Each.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Operations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Tracer.hpp"

namespace utopia {

    template <class InnerExpr, class Operation, class Traits, int Backend>
    class Eval<Reduce<InnerExpr, Operation>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;
        using Expr = utopia::Reduce<InnerExpr, Operation>;

        template <class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Expr &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<InnerExpr, Traits>::apply(expr.expr()).reduce(expr.operation());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Dot<Left, Right>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;

        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Dot<Left, Right> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).dot(
            //         Eval<Left,  Traits>::apply(expr.expr().left()),
            //         Eval<Right, Traits>::apply(expr.expr().right())
            // );

            result = Eval<Left, Traits>::apply(expr.expr().left()).dot(Eval<Right, Traits>::apply(expr.expr().right()));

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Norm<InnerExpr, 2>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;

        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Norm<InnerExpr, 2> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<InnerExpr, Traits>::apply(expr.expr()).norm2();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Norm<InnerExpr, 1>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;

        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Norm<InnerExpr, 1> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).norm1(
            //         Eval<Expr, Traits>::apply(expr.expr())
            // );

            result = Eval<InnerExpr, Traits>::apply(expr.expr()).norm1();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Norm<InnerExpr, INFINITY_NORM_TAG>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;

        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Norm<InnerExpr, INFINITY_NORM_TAG> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).norm_infty(
            //         Eval<Expr, Traits>::apply(expr.expr())
            // );

            result = Eval<InnerExpr, Traits>::apply(expr.expr()).norm_infty();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<ReduceApproxEqual<Left, Right>, Traits, Backend> {
    public:
        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static bool apply(const ReduceApproxEqual<Left, Right> &expr) {
            bool result;
            UTOPIA_TRACE_BEGIN(expr);

            // result = UTOPIA_BACKEND(Traits).compare(
            //         Eval<Left,  Traits>::apply(expr.expr().left()),
            //         Eval<Right, Traits>::apply(expr.expr().right()),
            //         expr.expr().operation());

            result = Eval<Left, Traits>::apply(expr.expr().left())
                         .equals(Eval<Right, Traits>::apply(expr.expr().right()), expr.expr().operation().tol());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    // TODO
    // template<class Expr, class Traits, int Backend>
    // class Eval<Reduce< Diag<Expr>, Plus>, Traits, Backend> {
    // public:
    //     typedef typename Traits::Scalar Scalar;
    //     inline static Number<Scalar> apply(const Reduce< Diag<Expr>, Plus> &expr)
    //     {
    //         Scalar result;
    //         UTOPIA_TRACE_BEGIN(expr);

    //         result = UTOPIA_BACKEND(Traits).trace(
    //                 Eval<Expr, Traits>::apply(expr.expr().expr())
    //         );

    //         UTOPIA_TRACE_END(expr);
    //         return result;
    //     }
    // };

    template <class InnerExpr, class Traits, int Backend>
    class Eval<Trace<InnerExpr>, Traits, Backend> {
    public:
        using Scalar = typename Traits::Scalar;

        // FIXME this lazy and not safe
        template <class Expr, class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static Number<Scalar> apply(const Trace<InnerExpr> &expr) {
            Scalar result;
            UTOPIA_TRACE_BEGIN(expr);

            result = Eval<InnerExpr, Traits>::apply(expr.expr()).trace();

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template <class Derived, typename T, class Traits, int Backend>
    class Eval<Reduce<Tensor<Derived, 2>, PlusIsNonZero<T>>, Traits, Backend> {
    public:
        using Expr = utopia::Reduce<Tensor<Derived, 2>, PlusIsNonZero<T>>;
        using Vector = UTOPIA_VECTOR(Derived);

        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;

        template <class Result>
        inline static void apply(const Expr &expr, Result &result) {
            result = apply(expr);
        }

        inline static SizeType apply(const Expr &expr) {
            SizeType result = 0;
            UTOPIA_TRACE_BEGIN(expr);

            const Scalar tol = expr.operation().is_non_zero().tol();

            const auto &mat = expr.expr().derived();
            mat.map_reduce(
                // map
                UTOPIA_LAMBDA(const Scalar &x)->SizeType { return static_cast<SizeType>(x > tol); },
                // reduce
                UTOPIA_LAMBDA(const SizeType &left, const SizeType &right)->SizeType { return left + right; },
                // mpi operation
                Plus(),
                result);

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_REDUCE_HPP
