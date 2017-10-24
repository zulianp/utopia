#ifndef UTOPIA_FUNCTION_SPACE_HPP
#define UTOPIA_FUNCTION_SPACE_HPP 

#include "utopia_Base.hpp"

namespace utopia {

	template<class Derived>
	class FunctionSpace {
	public:
		DERIVED_CRT(Derived);
		CONST_DERIVED_CRT(Derived);
	};
}

#endif //UTOPIA_FUNCTION_SPACE_HPP
