#ifndef UTOPIA_FE_EVAL_INTERPOLATE_HPP
#define UTOPIA_FE_EVAL_INTERPOLATE_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {
	template<class Coefficient, class Fun, class Traits, int Backend>
	class FEEval<Interpolate<Coefficient, Fun>, Traits, Backend> {
	public:
		typedef utopia::Interpolate<TrialFunction<Fun> > Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::interpolate(expr, ctx) )
	    {
	    	return FEBackend<Backend>::interpolate(expr, ctx);
	    } 
	};

	template<class Coefficient, class Fun, class AssemblyContext>
	class FunctionalTraits<Interpolate<Coefficient, Fun>, AssemblyContext> {
	public:
		inline static int type(const Interpolate<Coefficient, Fun> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Fun, AssemblyContext>::type(expr.expr(), ctx);
		}

		inline static int order(const Interpolate<Fun> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Fun, AssemblyContext>::order(expr.expr(), ctx);
		}
	};
}

#endif //UTOPIA_FE_EVAL_INTERPOLATE_HPP
