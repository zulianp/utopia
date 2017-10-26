#ifndef UTOPIA_FE_EVAL_EMPTY_HPP
#define UTOPIA_FE_EVAL_EMPTY_HPP 

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Eval.hpp"

namespace utopia {
	template<class Expr, class Traits, int Backend>
	class FEEval : public Eval<Expr, Traits, Backend> { 
	public:
		//default fallback on eval if does not exists
		inline static auto apply(const Expr &expr, const AssemblyContext<Backend> &) -> decltype( Eval<Expr, Traits, Backend>::apply(expr) )
		{
			return Eval<Expr>::apply(expr);
		}
	};

	template<class Tensor, int Order, class Traits, int Backend>
	class FEEval< Wrapper<Tensor, Order>, Traits, Backend> { 
	public:
		//default fallback on eval if does not exists
		inline static const Wrapper<Tensor, Order> &apply(const Wrapper<Tensor, Order> &expr, const AssemblyContext<Backend> &)
		{
			return expr;
		}
	};
}

#endif //UTOPIA_FE_EVAL_EMPTY_HPP
