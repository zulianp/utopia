//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
#define UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Right>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply<Left, Right> >::Type Result;

        inline static Result apply(const Multiply<Left, Right> &expr) {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).apply(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right()),
                    Multiplies(),
                    result
            );

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Multiply<Transposed<Left>, Right>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply < Transposed < Left>, Right> >::Type Result;

        inline static Result apply(const Multiply <Transposed<Left>, Right> &expr) {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).gemm(
                    1.0,
                    Eval<Left, Traits>::apply(expr.left().expr()),
                    Eval<Right, Traits>::apply(expr.right()),
                    true,
                    false,
                    0.0,
                    result);

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Transposed<Right> >, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply<Left, Transposed<Right> > >::Type Result;

        inline static Result apply(const Multiply<Left, Transposed<Right> > &expr)         {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).gemm(
                    1.0,
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr()),
                    false,
                    true,
                    0.0,
                    result);

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Transposed<Left>, Transposed<Right> >, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply<Transposed<Left>, Transposed<Right> > >::Type Result;

        inline static Result apply(const Multiply< Transposed<Left>, Transposed<Right> > &expr)
        {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).gemm(
                    1.0,
                    Eval<Left,  Traits>::apply(expr.left().expr()),
                    Eval<Right, Traits>::apply(expr.right().expr()),
                    true,
                    true,
                    0.0,
                    result);

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Transposed< Multiply<Left, Right> >, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Transposed< Multiply<Left, Right> > >::Type Result;

        inline static Result apply(const Transposed< Multiply<Left, Right> > &expr)
        {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            bool ok = UTOPIA_BACKEND(Traits).gemm(
                    1.0,
                    Eval<Right, Traits>::apply(expr.expr().right()),
                    Eval<Left, Traits>::apply(expr.expr().left()),
                    true,
                    true,
                    0.0,
                    result);

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
