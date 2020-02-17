#ifndef UTOPIA_PETSC_EVAL_FACTORY_HPP
#define UTOPIA_PETSC_EVAL_FACTORY_HPP

namespace utopia {
    template<class Index, class Traits>
    class Eval< Ghosts<Index>, Traits, PETSC> {
    public:
        typedef typename TypeAndFill<Traits, Ghosts<Index> >::Type Return;

        inline static Return apply(const Ghosts<Index> &expr) {
            Return ret;

            UTOPIA_TRACE_BEGIN(expr);

            ret.ghosted(
                expr.local_size(), expr.global_size(), expr.index()
            );

            UTOPIA_TRACE_END(expr);
            return ret;
        }
    };

    // template<class Left, class Index, class Traits>
    // class Eval< Construct<Left, Ghosts<Index> >, Traits, PETSC> {
    // public:
    //     typedef utopia::Construct<Left, Ghosts<Index> > Expr;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &r = expr.right();

    //         Eval<Left, Traits, PETSC>::apply(expr.left()).ghosted(
    //             r.local_size(), r.global_size(), r.index()
    //         );

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };


    template<class Left, class Index, class Traits>
    class Eval< Assign<Left, Ghosts<Index> >, Traits, PETSC> {
    public:
        typedef utopia::Assign<Left, Ghosts<Index> > Expr;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &r = expr.right();

            Eval<Left, Traits, PETSC>::apply(expr.left()).ghosted(
                r.local_size(), r.global_size(), r.index()
            );

            UTOPIA_TRACE_END(expr);
        }
    };
}

#endif //UTOPIA_PETSC_EVAL_FACTORY_HPP
