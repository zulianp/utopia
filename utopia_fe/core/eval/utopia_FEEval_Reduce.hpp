#ifndef UTOPIA_FE_EVAL_REDUCE_HPP
#define UTOPIA_FE_EVAL_REDUCE_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {

	template<class Expr, class AssemblyContext>
	class FunctionalTraits<Reduce<Expr, Plus>, AssemblyContext> {
	public:
		inline static int type(const Reduce<Expr, Plus> &expr, const AssemblyContext &ctx)  
		{ 
			return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);
		}

		inline static int order(const Reduce<Expr, Plus> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx);
		}
	};

}

#endif //UTOPIA_FE_EVAL_REDUCE_HPP
