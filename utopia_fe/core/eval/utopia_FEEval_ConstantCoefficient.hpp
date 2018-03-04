#ifndef UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP
#define UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP

namespace utopia {
	template<typename T, int Order, class Traits, int Backend, int IsQuadData>
	class FEEval<ConstantCoefficient<T, Order>, Traits, Backend, IsQuadData> { 
	public:
		typedef utopia::ConstantCoefficient<T, Order> Expr;

		inline static const Expr &apply(const Expr &expr, const AssemblyContext<Backend> &) {
			return expr;
		}
	};

	template<class T, int Order, class AssemblyContext>
	class FunctionalTraits<ConstantCoefficient<T, Order>, AssemblyContext> {
	public:
		inline static int type(const ConstantCoefficient<T, Order> &expr, const AssemblyContext &)  
		{ 
			return CONSTANT_FUNCTION;
		}

		inline static int order(const ConstantCoefficient<T, Order> &expr, const AssemblyContext &) 
		{
			return 0;
		}
	};
}

#endif //UTOPIA_FE_EVAL_CONSTANT_COEFFICIENT_HPP
