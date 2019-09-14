//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
#define UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Eval_Binary.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Right>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply<Left, Right> >::Type Result;

        inline static Result apply(const Multiply<Left, Right> &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            EvalBinaryAux<Result>::apply(
                Eval<Left,  Traits>::apply(expr.left()),
                Eval<Right, Traits>::apply(expr.right()),
                Multiplies(),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Multiply<Transposed<Left>, Right>, Traits, Backend> {
    public:
        typedef utopia::Multiply< Transposed<Left>, Right>  Expr;

        typedef typename TypeAndFill<Traits, Expr>::Type Result;

        inline static Result apply(const Expr &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            auto &&l = Eval<Left,  Traits>::apply(expr.left().expr());
            auto &&r = Eval<Right, Traits>::apply(expr.right());
            
            l.transpose_multiply(
                r, 
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Transposed<Right> >, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Multiply<Left, Transposed<Right> > >::Type Result;

        inline static Result apply(const Multiply<Left, Transposed<Right> > &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).multiply(
            //     result,
            //     false,
            //     Eval<Left, Traits>::apply(expr.left()),
            //     true,
            //     Eval<Right, Traits>::apply(expr.right().expr())
            //     );


            Eval<Left, Traits>::apply(expr.left()).multiply_transpose(
                Eval<Right, Traits>::apply(expr.right().expr()),
                result
            );

            UTOPIA_TRACE_END(expr);
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

            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).multiply(
            //     result,
            //     true,
            //     Eval<Left,  Traits>::apply(expr.left().expr()),
            //     true,
            //     Eval<Right, Traits>::apply(expr.right().expr())
            //     );

            Eval<Left, Traits>::apply(expr.left().expr()).multiply(
                true,
                true,
                Eval<Right, Traits>::apply(expr.right().expr()),
                result
            );

            UTOPIA_TRACE_END(expr);
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

            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).multiply(
            //     result,
            //     true,
            //     Eval<Right, Traits>::apply(expr.expr().right()),
            //     true,
            //     Eval<Left, Traits>::apply(expr.expr().left())
            //     );

            Eval<Right, Traits>::apply(expr.expr().right()).multiply(
                true,
                true,
                Eval<Left, Traits>::apply(expr.expr().left()),
                result
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };


    // template<class Left, class Right, class Traits, int Backend>
    // class Eval<Multiply<std::vector<Left>, Right>, Traits, Backend> {
    // public:
    //     typedef typename TypeAndFill<Traits, Multiply<Left, Right> >::Type ResultElement;

    //     inline static std::vector<ResultElement> apply(const Multiply<std::vector<Left>, Right> &expr)
    //     {
    //         std::vector<ResultElement> result(expr.left().size());
    //         for(std::size_t i = 0; i < result.size(); ++i) {
    //             result[i] = Eval<Multiply<Left, Right>>::apply(expr.left()[i], expr.right());
    //         }

    //         return result;
    //     }
    // };

}

#endif //UTOPIA_UTOPIA_EVAL_MULTIPLY_HPP
