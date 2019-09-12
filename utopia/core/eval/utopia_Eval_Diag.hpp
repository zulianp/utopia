//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_DIAG_HPP
#define UTOPIA_UTOPIA_EVAL_DIAG_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply< Left, Diag<Right> >, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Left) apply(const Multiply<Left, Diag<Right> > &expr) {
            static_assert(Right::Order == 1, "Right has to be a vector");
            EXPR_TYPE(Traits, Left) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag_scale_right(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr()),
                    result);

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply< Diag<Left>, Right>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Right) apply(const Multiply< Diag<Left>, Right> &expr)
        {
            EXPR_TYPE(Traits, Right) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag_scale_left(
                    result,
                    Eval<Left,  Traits>::apply(expr.left().expr()),
                    Eval<Right, Traits>::apply(expr.right()));

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Assign /////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Diag<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Diag<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 1>, Diag<Right> >, Traits, Backend> {
    public:
        using LeftTensor = utopia::Tensor<Left, 1>;
        inline static bool apply(const Assign<LeftTensor, Diag<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).build_diag(
                Eval<LeftTensor, Traits>::apply(expr.left())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 2>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Tensor<Left, 2>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).diag(
            //         Eval<Tensor<Left, 2>,  Traits>::apply(expr.left()),
            //         Eval<Right, Traits>::apply(expr.right().expr().expr())
            // );

            Eval<Tensor<Left, 2>,  Traits>::apply(expr.left()).diag(
                        Eval<Right, Traits>::apply(expr.right().expr().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 1>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Tensor<Left, 1>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Left,  Traits>::apply(expr.right().expr().expr()));


            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Construct /////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Left, Diag<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Diag<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Tensor<Left, 2>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        typedef utopia::Tensor<Left, 2> WLeft;

        inline static bool apply(const Construct<WLeft, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<WLeft,  Traits>::apply(expr.left()),
                    Eval<Diag<Right>, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Tensor<Left, 1>, Diag<Right>>, Traits, Backend> {
    public:
        using LeftExpr  = utopia::Tensor<Left, 1>;
        using RightExpr = utopia::Diag<Right>;

        inline static bool apply(const Construct<LeftExpr, RightExpr> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).build_diag(
                Eval<LeftExpr,  Traits>::apply(expr.left())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Tensor<Left, 1>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Tensor<Left, 1>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Left,  Traits>::apply(expr.right().expr().expr()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class T, class Traits, int Backend>
    class Eval< Diag<Tensor<T, 2> >, Traits, Backend> {
    public:
        typedef utopia::Tensor<T, 2> WT;
        typedef utopia::Diag<WT> Expr;
        typedef typename Traits::Vector Result;

        inline static Result apply(const Expr &expr)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    result,
                    Eval<WT,  Traits>::apply(expr.expr())
                    );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_DIAG_HPP
