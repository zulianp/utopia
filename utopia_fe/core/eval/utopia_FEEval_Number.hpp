#ifndef UTOPIA_FE_EVAL_NUMBER_HPP
#define UTOPIA_FE_EVAL_NUMBER_HPP 

#include "utopia_FunctionalTraits.hpp"

namespace utopia {
	template<class T, class AssemblyContext>
	class FunctionalTraits<Number<T>, AssemblyContext> {
	public:
		inline static int type(const Number<T> &expr, const AssemblyContext &)  
		{ 
			return CONSTANT_FUNCTION;
		}

		inline static int order(const Number<T> &expr, const AssemblyContext &) 
		{
			return 0;
		}
	};
}

#endif //UTOPIA_FE_EVAL_NUMBER_HPP
