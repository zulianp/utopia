#ifndef UTOPIA_FE_FUNCTIONAL_TRAITS_HPP
#define UTOPIA_FE_FUNCTIONAL_TRAITS_HPP 

#include <iostream>

namespace utopia {
	static const int POLYNOMIAL_FUNCTION = 0;
	static const int EXPONENTIAL_FUNCTION = 1;
	static const int UNDEFINED_FUNCTION = 2;

	template<class Expr, class AssemblyContext>
	class FunctionalTraits {
	public:
#ifdef NDEBUG
		inline constexpr static int type(const Expr &,  const AssemblyContext &) { return POLYNOMIAL_FUNCTION; }
		inline constexpr static int order(const Expr &, const AssemblyContext &) { return 0; }
#else
		inline static int type(const Expr &expr,  const AssemblyContext &)
		{
			std::cout << "Unhandled FunctionalTraits type for " <<  expr.getClass() << std::endl;
			return POLYNOMIAL_FUNCTION;
		}

		inline static int order(const Expr &expr, const AssemblyContext &)
		{
			std::cout << "Unhandled FunctionalTraits order for " <<  expr.getClass() << std::endl;
			return 0; 
		}
#endif
		
	};

	template<class Derived, class AssemblyContext>
	inline int functional_order(const Expression<Derived> &expr, const AssemblyContext &ctx)
	{
		return FunctionalTraits<Derived, AssemblyContext>::order(expr.derived(), ctx);
	}

	template<class Derived, class AssemblyContext>
	inline int functional_type(const Expression<Derived> &expr, const AssemblyContext &ctx)
	{
		return FunctionalTraits<Derived, AssemblyContext>::type(expr.derived(), ctx);
	}
}

#endif //UTOPIA_FE_FUNCTIONAL_TRAITS_HPP
