#ifndef UTOPIA_FE_FORWARD_DECLARATIONS_HPP
#define UTOPIA_FE_FORWARD_DECLARATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

	class FEExpression;

	template<class> class Coefficient;
	template<class> class Curl;
	template<class> class DifferentialOperator;
	template<class> class Div;
	template<class> class FormExpressions;
	template<class> class FunctionSpace;
	template<class> class Grad;
	template<class> class Integral;
	template<class> class TestFunction;
	template<class> class TrialFunction;
	template<class> class BlockVar;

	template<class, class> class Interpolate;
	template<class, int> class TimeDerivative;
	// template<class...> class ProductFunctionSpace;
	template<class> class ProductFunctionSpace;
	
}

#endif //UTOPIA_FE_FORWARD_DECLARATIONS_HPP

