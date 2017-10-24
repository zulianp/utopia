#ifndef UTOPIA_TEST_FUNCTION_HPP
#define UTOPIA_TEST_FUNCTION_HPP 

#include "utopia_BasisFunction.hpp"
#include "utopia_FunctionSpace.hpp"

namespace utopia {
	template<class FunctionSpaceT>
	class TestFunction : public BasisFunction< TestFunction<FunctionSpaceT> > {
	public:
		static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;
	};

	template<class Derived>
	inline TestFunction<Derived> test(const FunctionSpace<Derived> &space)
	{
		//FIXME
		return TestFunction<Derived>();
	}	
}

#endif //UTOPIA_TEST_FUNCTION_HPP
