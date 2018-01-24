#ifndef UTOPIA_FE_EVAL_DIV_HPP
#define UTOPIA_FE_EVAL_DIV_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

	template<class Tensor, class Traits, int Backend>
	class FEEval<Divergence<Tensor>, Traits, Backend> {
	public:
		typedef utopia::Divergence<Tensor> Expr;
		typedef typename Traits::DivergenceType DivergenceT;

	    template<template<class> class Function, class Spaces>
	    inline static auto apply(
	    	const Divergence< Function<ProductFunctionSpace<Spaces> > > &expr,
	    	AssemblyContext<Backend> &ctx) -> DivergenceT
	    {
	    	const auto & space_ptr = expr.expr().space_ptr();
	    	return FEBackend<Backend>::div(expr.expr(), ctx);
	    }
	};

}

#endif //UTOPIA_FE_EVAL_DIV_HPP
