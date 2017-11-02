#ifndef UTOPIA_FE_EVAL_REDUCE_HPP
#define UTOPIA_FE_EVAL_REDUCE_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {

	template<class Traits, int Order>
	class InnerProductType {};

	template<class Traits>
	class InnerProductType<Traits, 0> {
	public:
		typedef typename Traits::Scalar Type;
	};

	template<class Traits>
	class InnerProductType<Traits, 1> {
	public:
		typedef typename Traits::Vector Type;
	};

	template<class Traits>
	class InnerProductType<Traits, 2> {
	public:
		typedef typename Traits::Matrix Type;
	};

	template<class Left, class Right, class Traits, int Backend>
	class FEEval< Reduce< Binary<Left, Right, EMultiplies>, Plus>, Traits, Backend> {
	public:
		typedef utopia::Reduce< Binary<Left, Right, EMultiplies>, Plus> Expr;
		static const int has_trial = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
		static const int has_test  = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
		
		typedef typename InnerProductType<Traits, has_trial + has_test>::Type Result;

		template<template<class> class Function, class Space>
	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> Result
	    {
	    	return FEBackend<Backend>::inner(
	    			FEEval<Left,  Traits, Backend>::apply(expr.expr().left(),  ctx),
	    			FEEval<Right, Traits, Backend>::apply(expr.expr().right(), ctx),
	    			ctx
	    		);
	    } 
	};

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
