#ifndef UTOPIA_FE_EVAL_INTERPOLATE_HPP
#define UTOPIA_FE_EVAL_INTERPOLATE_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {
	template<class Coefficient, class Fun, class Traits, int Backend, int IsQuadData>
	class FEEval<Interpolate<Coefficient, Fun>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Interpolate<Coefficient, Fun> Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::fun(expr, ctx) )
	    {
	    	return FEBackend<Backend>::fun(expr, ctx);
	    } 
	};

	template<class Coefficient, class Fun, class Traits, int Backend, int IsQuadData>
	class FEEval<Gradient<Interpolate<Coefficient, Fun>>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Gradient<Interpolate<Coefficient, Fun> > Expr;

	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::grad(expr.expr(), ctx) )
	    {
	    	return FEBackend<Backend>::grad(expr.expr(), ctx);
	    } 
	};

	template<class Coefficient, class Fun, class AssemblyContext>
	class FunctionalTraits<Interpolate<Coefficient, Fun>, AssemblyContext> {
	public:
		inline static int type(const Interpolate<Coefficient, Fun> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Fun, AssemblyContext>::type(expr.fun(), ctx);
		}

		inline static int order(const Interpolate<Coefficient, Fun> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Fun, AssemblyContext>::order(expr.fun(), ctx);
		}
	};
}

#endif //UTOPIA_FE_EVAL_INTERPOLATE_HPP
