#ifndef UTOPIA_FE_EVAL_MULTIPLY_HPP
#define UTOPIA_FE_EVAL_MULTIPLY_HPP 

#include "utopia_FEEval_Empty.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
	template<class Left, class Function,  class Traits, int Backend>
	class FEEval<Multiply<Left, Gradient<Function> >, Traits, Backend> {
	public:
		typedef Multiply<Left, Gradient<Function> > Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), ctx);
	    } 
	};
}

#endif //UTOPIA_FE_EVAL_MULTIPLY_HPP
