#ifndef UTOPIA_FE_EVAL_GRAD_HPP
#define UTOPIA_FE_EVAL_GRAD_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_fe_lang.hpp"

namespace utopia {

	template<class Tensor, class Traits, int Backend>
	class Eval<Grad<Tensor>, Traits, Backend> {
	public:
		typedef Grad<Tensor> Expr;
		typedef EXPR_TYPE(Traits, Expr) Result;

	    inline static Result apply(const Expr &expr) {
	    	return Result();
	    }
	};
}

#endif //UTOPIA_FE_EVAL_GRAD_HPP
