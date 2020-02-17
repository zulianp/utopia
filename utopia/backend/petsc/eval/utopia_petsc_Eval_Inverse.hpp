#ifndef UTOPIA_EVAL_INVERSE_PETSC_HPP
#define UTOPIA_EVAL_INVERSE_PETSC_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Inverse.hpp"

namespace utopia {

    // template<class Left, class Right, class Traits>
    // class Eval< Construct<Left, Inverse<Right> >, Traits, PETSC> {
    // public:
    //     typedef utopia::Construct<Left, Inverse<Right> > Expr;
    //     typedef typename TypeAndFill<Traits, Left>::Type Result;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);
    //         auto & left   = Eval<Left, Traits>::apply(expr.left());
    //         auto && right = Eval<Right, Traits>::apply(expr.right().expr());

    //         right.inverse(left);

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    template<class Left, class Right, class Traits>
    class Eval< Assign<Left, Inverse<Right> >, Traits, PETSC> {
    public:
        typedef utopia::Assign<Left, Inverse<Right> > Expr;
        typedef typename TypeAndFill<Traits, Left>::Type Result;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto & left   = Eval<Left, Traits>::apply(expr.left());
            auto && right = Eval<Right, Traits>::apply(expr.right().expr());

            right.inverse(left);

            UTOPIA_TRACE_END(expr);
        }
    };

}

#endif //UTOPIA_EVAL_INVERSE_PETSC_HPP
