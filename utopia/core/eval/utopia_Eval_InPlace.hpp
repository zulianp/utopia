//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_INPLACE_HPP
#define UTOPIA_UTOPIA_EVAL_INPLACE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<InPlace<Left, Right, Operation>, Traits, Backend> {
    public:

        inline static bool apply(const InPlace<Left, Right, Operation> &expr)
        {
            //FIXME connect to backend without tree transformation
            typedef utopia::Binary<Left, Right, Operation> TransformedExpr;
            typedef utopia::Assign<Left, TransformedExpr> AssignExpr;

            UTOPIA_TRACE_BEGIN(expr);

            bool out = Eval<AssignExpr, Traits>::apply(
                    AssignExpr(
                            expr.left(),
                            TransformedExpr(
                                    expr.left(),
                                    expr.right(),
                                    expr.operation()
                            )
                    )
            );

            UTOPIA_TRACE_END(expr);
            return out;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<InPlace<Left, Right, Plus>, Traits, Backend> {
    public:

        inline static bool apply(const InPlace<Left, Right, Plus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).axpy(
            //       Eval<Left, Traits>::apply(expr.left()),
            //       1.,
            //       Eval<Right, Traits>::apply(expr.right())
            //       );

            Eval<Left, Traits>::apply(expr.left()).axpy(1.0, Eval<Right, Traits>::apply(expr.right()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<InPlace<Left, Right, Minus>, Traits, Backend> {
    public:

        inline static bool apply(const InPlace<Left, Right, Minus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).axpy(
            //       Eval<Left, Traits>::apply(expr.left()),
            //       -1.,
            //       Eval<Right, Traits>::apply(expr.right())
            //       );


            
            Eval<Left, Traits>::apply(expr.left()).axpy(
                  -1.,
                  Eval<Right, Traits>::apply(expr.right())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, typename T, class Right, class Traits, int Backend>
    class Eval<InPlace<Left, Binary<Number<T>, Right, Multiplies>, Plus>, Traits, Backend> {
    public:
        typedef utopia::InPlace<Left, Binary<Number<T>, Right, Multiplies>, Plus> Expr;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).axpy(
            //       Eval<Left, Traits>::apply(expr.left()),
            //       static_cast<T>(expr.right().left()),
            //       Eval<Right, Traits>::apply(expr.right().right())
            //       );

            Eval<Left, Traits>::apply(expr.left()).axpy(expr.right().left(), Eval<Right, Traits>::apply(expr.right().right()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, typename T, class Right, class Traits, int Backend>
    class Eval<InPlace<Left, Binary<Number<T>, Right, Multiplies>, Minus>, Traits, Backend> {
    public:
        typedef utopia::InPlace<Left, Binary<Number<T>, Right, Multiplies>, Minus> Expr;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).axpy(
                  Eval<Left, Traits>::apply(expr.left()),
                  -static_cast<T>(expr.right().left()),
                  Eval<Right, Traits>::apply(expr.right().right())
                  );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<InPlace<Left, Number<Right>, Multiplies>, Traits, Backend> {
    public:
        inline static bool apply(const InPlace<Left, Number<Right>, Multiplies> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).scale(
            //     Eval<Left, Traits>::apply(expr.left()),
            //     static_cast<Right>(expr.right())
            // );

            Eval<Left, Traits>::apply(expr.left()).scale(expr.right());

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_INPLACE_HPP
