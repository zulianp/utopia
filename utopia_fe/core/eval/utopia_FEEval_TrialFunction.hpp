#ifndef UTOPIA_FE_EVAL_TRIAL_FUNCTION_HPP
#define UTOPIA_FE_EVAL_TRIAL_FUNCTION_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {

	template<class FunctionSpaceT, class Traits, int Backend, int IsQuadData>
	class FEEval<TrialFunction<FunctionSpaceT>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::TrialFunction<FunctionSpaceT> Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::fun(expr, ctx) )
	    {
	    	return FEBackend<Backend>::fun(expr, ctx);
	    } 
	};

}

#endif //UTOPIA_FE_EVAL_TRIAL_FUNCTION_HPP
