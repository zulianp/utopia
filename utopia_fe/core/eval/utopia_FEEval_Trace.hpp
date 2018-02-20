#ifndef UTOPIA_FE_EVAL_TRACE_HPP
#define UTOPIA_FE_EVAL_TRACE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

	template<class InnerExpr, class Traits, int Backend, int IsQuadData>
	class FEEval<Trace<InnerExpr>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Trace<InnerExpr> Expr;		

		inline static auto apply(
			const Expr &expr,
			AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::trace(FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx) )
		{
			return FEBackend<Backend>::trace(FEEval<InnerExpr, Traits, Backend, IsQuadData>::apply(expr.expr(), ctx), ctx);
		}  
	};

	template<class Inner, class AssemblyContext>
	class FunctionalTraits<Trace<Inner>, AssemblyContext> {
	public:
		inline static int type(const Trace<Inner> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Inner, AssemblyContext>::type(expr.expr(), ctx);
		}

		inline static int order(const Trace<Inner> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Inner, AssemblyContext>::order(expr.expr(), ctx);
		}
	};
}

#endif //UTOPIA_FE_EVAL_TRACE_HPP
