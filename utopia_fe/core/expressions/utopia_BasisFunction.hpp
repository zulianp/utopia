#ifndef UTOPIA_BASIS_FUNCTION_HPP
#define UTOPIA_BASIS_FUNCTION_HPP 

#include "utopia_Expression.hpp"

namespace utopia {
	template<class Derived>
	class BasisFunction : public Expression<Derived> {};
}

#endif //UTOPIA_BASIS_FUNCTION_HPP
