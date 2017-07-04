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

            UTOPIA_LOG_BEGIN(expr);

            auto & cleft       = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().implementation();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

        	ASSERT(&cleft != &right && "should never happen");
        	const bool ok = UTOPIA_BACKEND(Traits).apply(left, right, Multiplies(), cleft); ASSERT(ok);

			UTOPIA_LOG_END(expr);
        }
    };

    template<class CLeft, class Left, class Right, class Traits, int Backend>
    class Eval< Assign<CLeft, Multiply<Wrapper<Left, 2>, Right> >, Traits, Backend> {
    public:
    	typedef utopia::Assign<CLeft, Multiply<Wrapper<Left, 2>, Right> > Expr;
        typedef typename TypeAndFill<Traits, CLeft>::Type Result;

        inline static void apply(const Expr &expr) {

            UTOPIA_LOG_BEGIN(expr);

            auto & cleft       = Eval<CLeft, Traits>::apply(expr.left());
            const auto & left  = expr.right().left().implementation();
            auto && right      = Eval<Right, Traits>::apply(expr.right().right());

        	if(&cleft != &right) {
        		const bool ok = UTOPIA_BACKEND(Traits).apply(left, right, Multiplies(), cleft);  ASSERT(ok);
        	} else {
        		Result result;
        		const bool ok = UTOPIA_BACKEND(Traits).apply(left, right, Multiplies(), result);  ASSERT(ok);
        		cleft = result;
        	}

            UTOPIA_LOG_END(expr);
        }
    };
}

#endif //UTOPIA_EVAL_CONSTRUCT_MULTIPLY_HPP
