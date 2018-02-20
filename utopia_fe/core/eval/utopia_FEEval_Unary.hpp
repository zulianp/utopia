#ifndef UTOPIA_FE_EVAL_UNARY_HPP
#define UTOPIA_FE_EVAL_UNARY_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"

namespace utopia {
	template<class Inner, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval< Unary<Inner, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Unary<Inner, Op> Expr;
		
	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype(
	    		FEBackend<Backend>::apply_unary(FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), expr.operation(), ctx)
	    	)
	    {
	    	return FEBackend<Backend>::apply_unary(FEEval<Inner, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), expr.operation(), ctx);
	    }  
	};


	template<class Inner, class Op, class AssemblyContext>
	class FunctionalTraits<Unary<Inner, Op>, AssemblyContext> {
	public:
		inline static int type(const Unary<Inner, Op> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Inner, AssemblyContext>::type(expr.expr(), ctx);
		}

		inline static int order(const Unary<Inner, Op> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Inner, AssemblyContext>::order(expr.expr(), ctx);
		}
	};

	template<class Inner, class AssemblyContext>
	class FunctionalTraits<Negate<Inner>, AssemblyContext> {
	public:
		inline static int type(const Negate<Inner> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Inner, AssemblyContext>::type(expr.expr(), ctx);
		}

		inline static int order(const Negate<Inner> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Inner, AssemblyContext>::order(expr.expr(), ctx);
		}
	};
}

#endif //UTOPIA_FE_EVAL_UNARY_HPP
