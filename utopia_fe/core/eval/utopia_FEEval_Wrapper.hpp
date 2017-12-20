#ifndef UTOPIA_FE_EVAL_WRAPPER_HPP
#define UTOPIA_FE_EVAL_WRAPPER_HPP 

#include "utopia_FunctionalTraits.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
	template<class T, int Order, class AssemblyContext>
	class FunctionalTraits<Wrapper<T, Order>, AssemblyContext> {
	public:
		inline static int type(const Wrapper<T, Order> &expr, const AssemblyContext &)  
		{ 
			return CONSTANT_FUNCTION;
		}

		inline static int order(const Wrapper<T, Order> &expr, const AssemblyContext &) 
		{
			return 0;
		}
	};
}

#endif //UTOPIA_FE_EVAL_WRAPPER_HPP
