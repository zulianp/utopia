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

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag_scale_right(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr()),
                    result);

            // [new backend map concept]
            // [fixme] Diag can be made a unary operator
            // [minimal] backend
            //(REMOVE)

            // [optimized] backend
            // UTOPIA_BACKEND(Traits).apply(result, left, Multiplies, DiagOp, right);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Multiply< Diag<Left>, Right>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Right) apply(const Multiply< Diag<Left>, Right> &expr)
        {
            EXPR_TYPE(Traits, Right) result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag_scale_left(
                    Eval<Left,  Traits>::apply(expr.left().expr()),
                    Eval<Right, Traits>::apply(expr.right()),
                    result);

            UTOPIA_LOG_END(expr);
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
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr())
            );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Wrapper<Left, 2>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Wrapper<Left, 2>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<Wrapper<Left, 2>,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr().expr())
            );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Wrapper<Left, 1>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Wrapper<Left, 1>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Left,  Traits>::apply(expr.right().expr().expr()));

            //FIXME error handling

            UTOPIA_LOG_END(expr);
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
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr())
            );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Wrapper<Left, 2>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        typedef utopia::Wrapper<Left, 2> WLeft;

        inline static bool apply(const Construct<WLeft, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    Eval<WLeft,  Traits>::apply(expr.left()),
                    // Eval<Right, Traits>::apply(expr.right().expr().expr())
                    Eval<Diag<Right>, Traits>::apply(expr.right().expr())
            );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Wrapper<Left, 1>, Diag< Diag<Right> > >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Wrapper<Left, 1>, Diag< Diag<Right> > > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Left,  Traits>::apply(expr.right().expr().expr()));

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Tensor, class Traits, int Backend>
    class Eval< Diag<Wrapper<Tensor, 2> >, Traits, Backend> {
    public:
        typedef utopia::Wrapper<Tensor, 2> WTensor;
        typedef utopia::Diag<WTensor> Expr;
        typedef typename Traits::Vector Result;

        inline static Result apply(const Expr &expr)
        {
            Result result;

            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).diag(
                    result,
                    Eval<WTensor,  Traits>::apply(expr.expr())
                    );

			UTOPIA_LOG_END(expr);
            return result;
        }
    };

    // template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Multiply< Left, Unary<Right, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static EXPR_TYPE(Traits, Left) apply(const Multiply<Left, Unary<Right, DiagOp> > &expr) {
    //             static_assert(Right::Order == 1, "Right has to be a vector");
    //             EXPR_TYPE(Traits, Left) result;

    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag_scale_right(
    //                     Eval<Left,  Traits>::apply(expr.left()),
    //                     Eval<Right, Traits>::apply(expr.right().expr()),
    //                     result);

    //             // [new backend map concept]
    //             // [fixme] Unary can be made a unary operator
    //             // [minimal] backend
    //             //(REMOVE)

    //             // [optimized] backend
    //             // UTOPIA_BACKEND(Traits).apply(result, left, Multiplies, UnaryOp, right);

    //             UTOPIA_LOG_END(expr);
    //             return result;
    //         }
    //     };

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Multiply< Unary<Left, DiagOp>, Right>, Traits, Backend> {
    //     public:
    //         inline static EXPR_TYPE(Traits, Right) apply(const Multiply< Unary<Left, DiagOp>, Right> &expr)
    //         {
    //             EXPR_TYPE(Traits, Right) result;

    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag_scale_left(
    //                     Eval<Left,  Traits>::apply(expr.left().expr()),
    //                     Eval<Right, Traits>::apply(expr.right()),
    //                     result);

    //             UTOPIA_LOG_END(expr);
    //             return result;
    //         }
    //     };

    //     //////////////////////////////////////////////////////////////////////////////////////////
    //     ///////////////////////////////////////// Assign /////////////////////////////////////////
    //     //////////////////////////////////////////////////////////////////////////////////////////

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Assign<Left, Unary<Right, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static bool apply(const Assign<Left, Unary<Right, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag(
    //                     Eval<Left,  Traits>::apply(expr.left()),
    //                     Eval<Right, Traits>::apply(expr.right().expr())
    //             );

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Assign<Wrapper<Left, 2>, Unary< Unary<Right, DiagOp>, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static bool apply(const Assign<Wrapper<Left, 2>, Unary< Unary<Right, DiagOp>, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag(
    //                     Eval<Wrapper<Left, 2>,  Traits>::apply(expr.left()),
    //                     Eval<Right, Traits>::apply(expr.right().expr().expr())
    //             );

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Assign<Wrapper<Left, 1>, Unary< Unary<Right, DiagOp>, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static bool apply(const Assign<Wrapper<Left, 1>, Unary< Unary<Right, DiagOp>, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).assign(
    //                     Eval<Left,  Traits>::apply(expr.left()),
    //                     Eval<Left,  Traits>::apply(expr.right().expr().expr()));

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     //////////////////////////////////////////////////////////////////////////////////////////
    //     ///////////////////////////////////////// Construct /////////////////////////////////////////
    //     //////////////////////////////////////////////////////////////////////////////////////////

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Construct<Left, Unary<Right, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static bool apply(const Construct<Left, Unary<Right, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag(
    //                     Eval<Left,  Traits>::apply(expr.left()),
    //                     Eval<Right, Traits>::apply(expr.right().expr())
    //             );

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Construct<Wrapper<Left, 2>, Unary< Unary<Right, DiagOp>, DiagOp> >, Traits, Backend> {
    //     public:
    //         typedef utopia::Wrapper<Left, 2> WLeft;

    //         inline static bool apply(const Construct<WLeft, Unary< Unary<Right, DiagOp>, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag(
    //                     Eval<WLeft,  Traits>::apply(expr.left()),
    //                     // Eval<Right, Traits>::apply(expr.right().expr().expr())
    //                     Eval<Unary<Right, DiagOp>, Traits>::apply(expr.right().expr())
    //             );

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     template<class Left, class Right, class Traits, int Backend>
    //     class Eval<Construct<Wrapper<Left, 1>, Unary< Unary<Right, DiagOp>, DiagOp> >, Traits, Backend> {
    //     public:
    //         inline static bool apply(const Construct<Wrapper<Left, 1>, Unary< Unary<Right, DiagOp>, DiagOp> > &expr)
    //         {
    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).assign(
    //                     Eval<Left,  Traits>::apply(expr.left()),
    //                     Eval<Left,  Traits>::apply(expr.right().expr().expr()));

    //             //FIXME error handling

    //             UTOPIA_LOG_END(expr);
    //             return true;
    //         }
    //     };

    //     template<class Tensor, class Traits, int Backend>
    //     class Eval< Unary<Wrapper<Tensor, 2>, DiagOp>, Traits, Backend> {
    //     public:
    //         typedef utopia::Wrapper<Tensor, 2> WTensor;
    //         typedef utopia::Unary<WTensor, DiagOp> Expr;
    //         typedef typename Traits::Vector Result;

    //         inline static Result apply(const Expr &expr)
    //         {
    //             Result result;

    //             UTOPIA_LOG_BEGIN(expr);

    //             UTOPIA_BACKEND(Traits).diag(
    //                     result,
    //                     Eval<WTensor,  Traits>::apply(expr.expr())
    //                     );

    //             UTOPIA_LOG_END(expr);
    //             return result;
    //         }
    //     };
}

#endif //UTOPIA_UTOPIA_EVAL_DIAG_HPP
