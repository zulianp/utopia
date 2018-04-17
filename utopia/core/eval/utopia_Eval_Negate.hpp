#ifndef UTOPIA_EVAL_NEGATE_HPP
#define UTOPIA_EVAL_NEGATE_HPP 

#include "utopia_Eval_Unary.hpp"

namespace utopia {
	template<class Expr, class Traits, int Backend>
	class Eval<Negate<Expr>, Traits, Backend> {
	public:
	  	typedef utopia::Unary<Expr, Minus> DelegateT;
	    typedef typename TypeAndFill<Traits, DelegateT>::Type Result;
	  
	    inline static Result apply(const Negate<Expr> &expr)
	    {
	        Result result;

	        UTOPIA_TRACE_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).apply_unary(
	            result,
	            expr.operation(),
	            Eval<DelegateT, Traits>::apply(expr.expr())
	            );

	        UTOPIA_TRACE_END(expr);
	        return result;
	    }
	};
}

#endif //UTOPIA_EVAL_NEGATE_HPP
