#ifndef UTOPIA_FE_EVAL_INTEGRAL_HPP
#define UTOPIA_FE_EVAL_INTEGRAL_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {
	template<class Inner, class Traits, int Backend>
	class FEEval< Integral<Inner>, Traits, Backend> {
	public:
		typedef utopia::Integral<Inner> Expr;
		
	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype(
	    		FEEval<Inner, Traits, Backend>::apply(expr.expr(), ctx)
	    	)
	    {
	    	return FEEval<Inner, Traits, Backend>::apply(expr.expr(), ctx);
	    }  
	};
}

#endif //UTOPIA_FE_EVAL_INTEGRAL_HPP
