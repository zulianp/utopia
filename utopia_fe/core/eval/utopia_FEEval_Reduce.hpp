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
		typedef typename Traits::Scalar Scalar;
	};

	template<class Traits>
	class InnerProductType<Traits, 1> {
	public:
		typedef typename Traits::Vector Vector;
	};

	template<class Traits>
	class InnerProductType<Traits, 2> {
	public:
		typedef typename Traits::Matrix Matrix;
	};

	// template<class Left, class Right, class Traits, int Backend>
	// class FEEval< Reduce< Binary<Left, Right, EMultiplies>, Plus> >, Traits, Backend> {
	// public:
	// 	typedef utopia::Reduce< Binary<Left, Right, EMultiplies>, Plus> > Expr;
		
	// 	typedef typename InnerProductType<Traits, I

	// 	template<template<class> class Function, class Space>
	//     inline static auto apply(
	//     	const Gradient< Function<Space> > &expr,
	//     	AssemblyContext<Backend> &ctx) -> GradientT
	//     {
	//     	return FEBackend<Backend>::grad(expr.expr(), ctx);
	//     } 

	//     class SubspaceVisitor {
	//     public:
	//     	template<class Subspace>
	//     	void operator()(const int i, const Subspace &space) 
	//     	{
	//     		std::cout << "visiting subspace: " << i << std::endl;
	//     	}
	//     };

	//     template<template<class> class Function, class Spaces>
	//     inline static auto apply(
	//     	const Gradient< Function<ProductFunctionSpace<Spaces> > > &expr,
	//     	AssemblyContext<Backend> &ctx) -> JacobianT
	//     {
	//     	const auto & space_ptr = expr.expr().space_ptr();
	//     	return FEBackend<Backend>::grad(expr.expr(), ctx);
	//     }

	// };



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
