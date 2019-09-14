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
    class Eval< Construct<Left, Transposed <Tensor<Right, 2> > >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, Transposed <Tensor<Right, 2> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto && left  = Eval<Left,  Traits>::apply(expr.left());
            auto && right = Eval<Tensor<Right, 2>, Traits>::apply(expr.right().expr());

            right.transpose(left);

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP
