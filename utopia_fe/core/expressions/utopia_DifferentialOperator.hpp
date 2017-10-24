#ifndef UTOPIA_FE_DIFFERENTIAL_OPERATOR_HPP
#define UTOPIA_FE_DIFFERENTIAL_OPERATOR_HPP 

#include "utopia_Expression.hpp"

namespace utopia {

	template<class Derived>
	class DifferentialOperator : public Expression<Derived> {
	public:
		/*virtual*/ ~DifferentialOperator() {}
	};

}


#endif //UTOPIA_FE_DIFFERENTIAL_OPERATOR_HPP
