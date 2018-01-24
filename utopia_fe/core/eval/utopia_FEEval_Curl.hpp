#ifndef UTOPIA_FE_EVAL_CURL_HPP
#define UTOPIA_FE_EVAL_CURL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"

namespace utopia {

	template<class Tensor, class Traits, int Backend>
	class FEEval<Curl<Tensor>, Traits, Backend> {
	public:
		typedef utopia::Curl<Tensor> Expr;
		typedef typename Traits::CurlType CurlT;

	    template<template<class> class Function, class Spaces>
	    inline static auto apply(
	    	const Curl< Function<ProductFunctionSpace<Spaces> > > &expr,
	    	AssemblyContext<Backend> &ctx) -> CurlT
	    {
	    	const auto & space_ptr = expr.expr().space_ptr();
	    	return FEBackend<Backend>::curl(expr.expr(), ctx);
	    }
	};

}

#endif //UTOPIA_FE_EVAL_CURL_HPP
