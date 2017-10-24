#ifndef UTOPIA_FE_EVAL_GRAD_HPP
#define UTOPIA_FE_EVAL_GRAD_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

	template<class Tensor, class Traits, int Backend>
	class FEEval<Gradient<Tensor>, Traits, Backend> {
	public:
		typedef utopia::Gradient<Tensor> Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::grad(expr.expr(), ctx) )
	    {
	    	return FEBackend<Backend>::grad(expr.expr(), ctx);
	    } 
	};
}

#endif //UTOPIA_FE_EVAL_GRAD_HPP
