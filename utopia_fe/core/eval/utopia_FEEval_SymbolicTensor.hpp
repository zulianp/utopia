#ifndef UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP
#define UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP 

#include "utopia_Factory.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"

namespace utopia {

	template<class Type, int Order, class AssemblyContext>
	class FunctionalTraits<SymbolicTensor<Type, Order>, AssemblyContext> {
	public:
		inline static int type(const SymbolicTensor<Type, Order> &expr, const AssemblyContext &ctx)  
		{ 
			return CONSTANT_FUNCTION;
		}

		inline static int order(const SymbolicTensor<Type, Order> &expr, const AssemblyContext &ctx) 
		{
			return 0;
		}
	};
}

#endif //UTOPIA_FE_EVAL_SYMBOLIC_TENSOR_HPP
