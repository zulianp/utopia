//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_FACTORY_HPP
#define UTOPIA_UTOPIA_EVAL_FACTORY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval< Assign<View<Left>, Factory<Right, Order> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<View<Left>, Factory<Right, Order> > &expr)
        {
            const auto &left = expr.left();
            UTOPIA_BACKEND(Traits).assignToRange(
                    Eval<Left, Traits>::apply(left.expr()),
                    expr.right().type(),
                    row_range(left),
                    col_range(left)
            );

            //FIXME error handling
            return true;
        }
    };

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval< Construct<View<Left>, Factory<Right, Order> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<View<Left>, Factory<Right, Order> > &expr)
        {
            const auto &left = expr.left();
            UTOPIA_BACKEND(Traits).assignToRange(
                    Eval<Left, Traits>::apply(left.expr()),
                    expr.right().type(),
                    row_range(left),
                    col_range(left)
            );

            //FIXME error handling
            return true;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Right, Order> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Factory<Right, Order> > &expr)
        {
            UTOPIA_BACKEND(Traits).build(
                    Eval<Left, Traits>::apply(expr.left()),
                    expr.right().size(),
                    expr.right().type()
            );

            //FIXME error handling
            return true;
        }
    };

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<Right, Order> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Factory<Right, Order> > &expr)
        {
            UTOPIA_BACKEND(Traits).build(
                    Eval<Left, Traits>::apply(expr.left()),
                    expr.right().size(),
                    expr.right().type()
            );

            //FIXME error handling
            return true;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Type, int Order, class Traits, int Backend>
    class Eval< Factory<Type, Order>, Traits, Backend> {
    public:
        inline static typename TypeAndFill<Traits, Factory<Type, Order> >::Type apply(const Factory<Type, Order> &expr) {
            typename TypeAndFill<Traits, Factory<Type, Order> >::Type ret;
            UTOPIA_BACKEND(Traits).build(ret, expr.size(), expr.type());
            return ret;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_FACTORY_HPP
