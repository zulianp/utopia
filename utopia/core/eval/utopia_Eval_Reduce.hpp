//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_REDUCE_HPP
#define UTOPIA_UTOPIA_EVAL_REDUCE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template <class Expr, class Operation, class Traits, int Backend>
    class Eval<Reduce<Expr, Operation>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Reduce<Expr, Operation> &expr)
        {
            Scalar result;
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).reduce(
                    Eval<Expr, Traits>::apply(expr.expr()),
                    expr.operation()
            );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Reduce<Binary<Left, Right, EMultiplies>, Plus>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr)
        {
            Scalar result;
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).dot(
                    Eval<Left,  Traits>::apply(expr.expr().left()),
                    Eval<Right, Traits>::apply(expr.expr().right())
            );

            UTOPIA_LOG_END(expr);
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
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).norm2(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Norm<Expr, 1>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Norm<Expr, 1> &expr) {
            Scalar result;
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).norm1(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Expr, class Traits, int Backend>
    class Eval<Norm<Expr, INFINITY_NORM_TAG>, Traits, Backend> {
    public:
        typedef typename Traits::Scalar Scalar;

        inline static Scalar apply(const Norm<Expr, INFINITY_NORM_TAG> &expr) {
            Scalar result;
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).norm_infty(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Reduce<Binary<Left, Right, ApproxEqual>, And>, Traits, Backend> {
    public:
        inline static bool apply(const Reduce<Binary<Left, Right, ApproxEqual>, And> &expr) {
            bool result;
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).compare(
                    Eval<Left,  Traits>::apply(expr.expr().left()),
                    Eval<Right, Traits>::apply(expr.expr().right()),
                    expr.expr().operation());

            UTOPIA_LOG_END(expr);
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
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).trace(
                    Eval<Expr, Traits>::apply(expr.expr().expr())
            );

            UTOPIA_LOG_END(expr);
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
            UTOPIA_LOG_BEGIN(expr);

            result = UTOPIA_BACKEND(Traits).trace(
                    Eval<Expr, Traits>::apply(expr.expr())
            );

            UTOPIA_LOG_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_REDUCE_HPP
