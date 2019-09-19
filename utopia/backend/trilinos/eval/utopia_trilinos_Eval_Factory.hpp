#ifndef UTOPIA_TRILINOS_EVAL_FACTORY_HPP
#define UTOPIA_TRILINOS_EVAL_FACTORY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    template<class Index, class Traits>
    class Eval< Ghosts<Index>, Traits, TRILINOS> {
    public:
        typedef typename TypeAndFill<Traits, Ghosts<Index> >::Type Return;

        inline static Return apply(const Ghosts<Index> &expr) {
            Return ret;

            UTOPIA_TRACE_BEGIN(expr);

            ret.ghosted(expr.local_size(), expr.global_size(), expr.index());

            UTOPIA_TRACE_END(expr);
            return ret;
        }
    };

    template<class Left, class Index, class Traits>
    class Eval< Construct<Left, Ghosts<Index> >, Traits, TRILINOS> {
    public:
        typedef utopia::Construct<Left, Ghosts<Index> > Expr;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            auto &&g = expr.right();

            Eval<Left, Traits, TRILINOS>::apply(expr.left()).ghosted(g.local_size(), g.global_size(), g.index());

            UTOPIA_TRACE_END(expr);
        }
    };

    template<class Left, class Index, class Traits>
    class Eval< Assign<Left, Ghosts<Index> >, Traits, TRILINOS> {
    public:
        typedef utopia::Assign<Left, Ghosts<Index> > Expr;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).build_ghosts(
                expr.right().local_size(),
                expr.right().global_size(),
                expr.right().index(),
                Eval<Left, Traits, TRILINOS>::apply(expr.left())
            );

            UTOPIA_TRACE_END(expr);
        }
    };
}

#endif //UTOPIA_TRILINOS_EVAL_FACTORY_HPP
