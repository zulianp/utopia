#ifndef UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP
#define UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP

namespace utopia {
	template<typename T, int Order, class Traits, int Backend>
	class FEEval<ConstantCoefficient<T, Order>, Traits, Backend> { 
	public:
		typedef utopia::ConstantCoefficient<T, Order> Expr;

		inline static const Expr &apply(const Expr &expr, const AssemblyContext<Backend> &) {
			return expr;
		}
	};
}

#endif //UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP
