#ifndef UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP
#define UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {
	template<class Out, class Fun, class Traits, int Backend>
	class FEEval<ContextFunction<Out, Fun>, Traits, Backend> {
	public:
		typedef utopia::ContextFunction<Out, Fun> Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> Out
	    {
	    	return expr.eval(ctx);
	    } 
	};
}


#endif //UTOPIA_FE_EVAL_CONTEXT_FUNCTION_HPP
