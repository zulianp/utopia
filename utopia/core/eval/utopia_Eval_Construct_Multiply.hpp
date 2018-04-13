#ifndef UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP
#define UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP


#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class CLeft, class Left, class Right, class Traits, int Backend>
    class Eval< Construct<CLeft, Multiply<Wrapper<Left, 2>, Right> >, Traits, Backend> {
    public:
    	typedef utopia::Construct<CLeft, Multiply<Wrapper<Left, 2>, Right> > Expr;
        typedef typename TypeAndFill<Traits, CLeft>::Type Result;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto & result      = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().implementation();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

        	ASSERT(&result != &right && "should never happen");

        	UTOPIA_BACKEND(Traits).apply_binary(result, left, Multiplies(), right);

			UTOPIA_TRACE_END(expr);
        }
    };

    template<class CLeft, class Left, class Right, class Traits, int Backend>
    class Eval< Assign<CLeft, Multiply<Wrapper<Left, 2>, Right> >, Traits, Backend> {
    public:
    	typedef utopia::Assign<CLeft, Multiply<Wrapper<Left, 2>, Right> > Expr;
        typedef typename TypeAndFill<Traits, CLeft>::Type Result;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto & result      = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().implementation();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

        	if(&result != &right) {
        		UTOPIA_BACKEND(Traits).apply_binary(result, left, Multiplies(), right);
        	} else {
        		typename std::remove_const<typename std::remove_reference<decltype(result)>::type>::type temp;
        		UTOPIA_BACKEND(Traits).apply_binary(temp, left, Multiplies(),right);
                result = temp;
        	}

            UTOPIA_TRACE_END(expr);
        }
    };
}

#endif //UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP
