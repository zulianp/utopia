#ifndef UTOPIA_EVAL_SELECT_HPP
#define UTOPIA_EVAL_SELECT_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    

    //TODO
    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Assign<Left, Select<Right, 1> >, Traits, Backend> {
    // public:
    //     inline static bool apply(const Assign<Left, Select<Right, 1> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).select(
    //             Eval<Left,  Traits>::apply(expr.left()),
    //             Eval<Right, Traits>::apply(expr.right().expr()),
    //             expr.right().index()
    //             );

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };


    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Construct<Left, Select<Right, 1> >, Traits, Backend> {
    // public:
    //     inline static bool apply(const Construct<Left, Select<Right, 1> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Right, Traits>::apply(expr.right().expr()).select(
    //             expr.right().index(),
    //             Eval<Left,  Traits>::apply(expr.left())
    //             );


    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };

    template<class Expr, class Traits, int Backend>
    class Eval< Select<Expr, 1>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Expr>::Type Result;

        inline static Result apply(const Select<Expr, 1> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            Result result;

            Eval<Expr, Traits>::apply(expr.expr()).select(
                expr.index(),
                result
                );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };


    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign<Left, Select<Right, 2> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Select<Right, 2> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).select(
                expr.right().row_index(),
                expr.right().col_index(),
                Eval<Left,  Traits>::apply(expr.left())
                );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Construct<Left, Select<Right, 2> >, Traits, Backend> {
    // public:
    //     inline static bool apply(const Construct<Left, Select<Right, 2> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Right, Traits>::apply(expr.right().expr()).select(
    //             expr.right().row_index(),
    //             expr.right().col_index(),
    //             Eval<Left,  Traits>::apply(expr.left())
    //             );

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };


    template<class Expr, class Traits, int Backend>
    class Eval< Select<Expr, 2>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Expr>::Type Result;

        inline static Result apply(const Select<Expr, 2> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            Result result;

            Eval<Expr, Traits>::apply(expr.expr()).select(
                expr.row_index(),
                expr.col_index(),
                result
                );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_EVAL_SELECT_HPP
