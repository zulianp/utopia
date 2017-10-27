#ifndef UTOPIA_FE_EVAL_BINARY_HPP
#define UTOPIA_FE_EVAL_BINARY_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {

	template<class Left, class Space, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, TrialFunction<Space>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, TrialFunction<Space>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};

	template<class Left, class Space, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, TestFunction<Space>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, TestFunction<Space>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};

	template<class Left, class Function, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, Gradient<Function>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, Gradient<Function>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};


	template<class Left, class Right, class Op, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, Op>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, Op> &expr, const AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, Op> &expr, const AssemblyContext &ctx) 
		{
			return std::max(FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx), 
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx));
		}
	};

	template<class Left, class Right, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, EMultiplies>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, EMultiplies> &expr, const AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, EMultiplies> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
		}
	};


	template<class Left, class Right, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, Multiplies>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
		}
	};
}


#endif //UTOPIA_FE_EVAL_BINARY_HPP
