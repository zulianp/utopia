#ifndef UTOPIA_BASIS_FUNCTION_HPP
#define UTOPIA_BASIS_FUNCTION_HPP 

#include "utopia_Expression.hpp"
#include "utopia_FEExpression.hpp"

#include <memory>

namespace utopia {
	template<class Derived, class FunctionSpaceT>
	class BasisFunction : public Expression<Derived>, public FEExpression {
	public:
		BasisFunction(const FunctionSpaceT &space)
		: space_ptr_(std::make_shared<FunctionSpaceT>(space)) 
		{}

		BasisFunction(std::shared_ptr<FunctionSpaceT> &space_ptr)
		: space_ptr_(space_ptr) 
		{}

		BasisFunction(FunctionSpaceT &&space)
		: space_ptr_(std::make_shared<FunctionSpaceT>(std::move(space)))
		{}

		inline std::shared_ptr<FunctionSpaceT> space_ptr() const
		{
			return space_ptr_;
		}

		std::shared_ptr<FunctionSpaceT> space_ptr_;
	};
}

#endif //UTOPIA_BASIS_FUNCTION_HPP
