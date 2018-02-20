#ifndef UTOPIA_FE_EVAL_GRAD_HPP
#define UTOPIA_FE_EVAL_GRAD_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"


namespace utopia {

	template<class Tensor, class Traits, int Backend, int IsQuadData>
	class FEEval<Gradient<Tensor>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Gradient<Tensor> Expr;
		typedef typename Traits::GradientType GradientT;
		typedef typename Traits::JacobianType JacobianT;

		template<template<class> class Function, class Space>
	    inline static auto apply(
	    	const Gradient< Function<Space> > &expr,
	    	AssemblyContext<Backend> &ctx) -> GradientT
	    {
	    	return FEBackend<Backend>::grad(expr.expr(), ctx);
	    } 

	    class SubspaceVisitor {
	    public:
	    	template<class Subspace>
	    	void operator()(const int i, const Subspace &space) 
	    	{
	    		std::cout << "visiting subspace: " << i << std::endl;
	    	}
	    };

	    template<template<class> class Function, class Spaces>
	    inline static auto apply(
	    	const Gradient< Function<ProductFunctionSpace<Spaces> > > &expr,
	    	AssemblyContext<Backend> &ctx) -> JacobianT
	    {
	    	const auto & space_ptr = expr.expr().space_ptr();
	    	return FEBackend<Backend>::grad(expr.expr(), ctx);
	    }

	};

}

#endif //UTOPIA_FE_EVAL_GRAD_HPP
