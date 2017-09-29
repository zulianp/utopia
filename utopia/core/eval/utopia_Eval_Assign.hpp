//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP

#include "utopia_Eval_Empty.hpp"


namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Right>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Right> &expr) {
            UTOPIA_LOG_BEGIN(expr);
            UTOPIA_BACKEND(Traits).assign(Eval<Left,  Traits>::apply(expr.left()),
                                          Eval<Right, Traits>::apply(expr.right()) );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign< View<Left>, Right>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<View<Left>, Right> &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            const auto &left = expr.left();
            auto rr = row_range(left);
            auto cr = col_range(left);

            UTOPIA_BACKEND(Traits).assign_to_range(Eval<Left,  Traits>::apply(expr.left().expr()),
                                                 Eval<Right, Traits>::apply(expr.right()),
                                                 rr, cr);
            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign< View< Wrapper<Left, 1> >, Right>, Traits, Backend> {
    public:
        typedef utopia::Wrapper<Left, 1> LeftWrapper;

        inline static bool apply(const Assign<View<LeftWrapper>, Right> &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            const auto &left = expr.left();
            auto rr = row_range(left);
            auto cr = col_range(left);

            UTOPIA_BACKEND(Traits).assign_to_range(Eval<LeftWrapper, Traits>::apply(expr.left().expr()),
                                                 Eval<Right, Traits>::apply(expr.right()),
                                                 rr, cr);

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign<Left, Transposed <Wrapper<Right, 2> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Transposed <Wrapper<Right, 2> > > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign_transposed(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Wrapper<Right, 2>, Traits>::apply(expr.right().expr())
            );

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign<Left, View<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, View<Right> > &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign_from_range(
                    Eval<Left,  Traits>::apply(expr.left()),
                    Eval<Right, Traits>::apply(expr.right().expr()),
                    row_range(expr.right()),
                    col_range(expr.right())
            );

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
