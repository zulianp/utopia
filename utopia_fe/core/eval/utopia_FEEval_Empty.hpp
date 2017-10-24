#ifndef UTOPIA_FE_EVAL_EMPTY_HPP
#define UTOPIA_FE_EVAL_EMPTY_HPP 

namespace utopia {
	template<class Expr, class Traits, int Backend>
	class FEEval : public Eval<Expr, Traits, Backend> { };
}

#endif //UTOPIA_FE_EVAL_EMPTY_HPP
