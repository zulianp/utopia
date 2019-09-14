#ifndef UTOPIA_PETSC_EVAL_PARALLEL_HPP
#define UTOPIA_PETSC_EVAL_PARALLEL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Tracer.hpp"

namespace utopia {

    template<class Left, class Right, class Traits>
    class Eval<Construct<Left, LocalDiagBlock<Right> >, Traits, PETSC> {
    public:

        inline static void apply(const Construct<Left, LocalDiagBlock<Right> > & expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).diagonal_block(
            	Eval<Left,  Traits>::apply(expr.left())
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    template<class Left, class Right, class Traits>
    class Eval<Assign<Left, LocalDiagBlock<Right> >, Traits, PETSC> {
    public:

        inline static void apply(const Assign<Left, LocalDiagBlock<Right> > & expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).diagonal_block(
            	Eval<Left,  Traits>::apply(expr.left())
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //FIXME todo
    // template<class Left, class Right, class Traits>
    // class Eval<LocalRedistribute<Left, Right>, Traits, PETSC> {
    // public:

    //     inline static EXPR_TYPE(Traits, Left) apply(const LocalRedistribute<Left, Right> &expr)
    //     {
    //         EXPR_TYPE(Traits, Left) result;

    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).build_local_redistribute(
    //             result,
    //             Eval<Left,  Traits>::apply(expr.left()),
    //             Eval<Right, Traits>::apply(expr.right())
    //             );


    //         UTOPIA_TRACE_END(expr);
    //         return result;
    //     }
    // };


}

#endif //UTOPIA_PETSC_EVAL_PARALLEL_HPP
