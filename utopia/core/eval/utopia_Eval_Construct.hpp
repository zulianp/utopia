//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Tracer.hpp"


namespace utopia {

    // [new backend map concept]
    // [minimal] construct and assign can be merged to the same
    // [optimized] backend (find out a way to provide user with better info about construction or assignment)


    template<class Left, class Right, class Traits, int Backend>
    class Eval<Construct<Left, Right>, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Right> &expr) {

            UTOPIA_TRACE_BEGIN(expr);

            expr.left().construct( Eval<Right, Traits>::apply(expr.right()) );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct< Number<Left>, Right>, Traits, Backend> {
    public:
        inline static bool apply(const Construct< Number<Left>, Right> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            expr.left().construct( Eval<Right, Traits>::apply(expr.right()) );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Left, Transposed <Wrapper<Right, 2> > >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Transposed <Wrapper<Right, 2> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).assign_transposed(
            //         Eval<Left,  Traits>::apply(expr.left()),
            //         Eval<Wrapper<Right, 2>, Traits>::apply(expr.right().expr())
            // );


            auto && left  = Eval<Left,  Traits>::apply(expr.left());
            auto && right = Eval<Wrapper<Right, 2>, Traits>::apply(expr.right().expr());

            right.transpose(left);

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };










    ///////TODOs

    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Construct< View<Left>, Right>, Traits, Backend> {
    // public:
    //     inline static bool apply(const Construct<View<Left>, Right> &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &left = expr.left();
    //         auto rr = row_range(left);
    //         auto cr = col_range(left);

    //         UTOPIA_BACKEND(Traits).assign_to_range(Eval<Left,  Traits>::apply(expr.left().expr()),
    //                                              Eval<Right, Traits>::apply(expr.right()),
    //                                              rr, cr);

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };

    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Construct< View< Wrapper<Left, 1> >, Right>, Traits, Backend> {
    // public:
    // typedef utopia::Wrapper<Left, 1> LeftWrapper;

    //     inline static bool apply(const Construct<View<LeftWrapper>, Right> &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &left = expr.left();
    //         auto rr = row_range(left);
    //         auto cr = col_range(left);

    //         UTOPIA_BACKEND(Traits).assign_to_range(Eval<LeftWrapper, Traits>::apply(expr.left().expr()),
    //                                              Eval<Right, Traits>::apply(expr.right()),
    //                                              rr, cr);

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };


    // template<class Left, class Right, class Traits, int Backend>
    // class Eval< Construct<Left, View<Right> >, Traits, Backend> {
    // public:
    //     inline static bool apply(const Construct<Left, View<Right> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).assign_from_range(
    //                 Eval<Left,  Traits>::apply(expr.left()),
    //                 Eval<Right, Traits>::apply(expr.right().expr()),
    //                 row_range(expr.right()),
    //                 col_range(expr.right())
    //         );

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };
}

#endif //UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP
