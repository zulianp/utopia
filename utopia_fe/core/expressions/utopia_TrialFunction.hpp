#ifndef UTOPIA_TRIAL_FUNCTION_HPP
#define UTOPIA_TRIAL_FUNCTION_HPP 

#include "utopia_BasisFunction.hpp"
#include "utopia_FunctionSpace.hpp"

namespace utopia {
	template<class FunctionSpaceT>
	class TrialFunction : public BasisFunction< TrialFunction<FunctionSpaceT> > {
	public:
		static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;
	};

	template<class Derived>
	TrialFunction<Derived> trial(const FunctionSpace<Derived> &space)
	{
		//FIXME
		return TrialFunction<Derived>();
	}
}

#endif //UTOPIA_TRIAL_FUNCTION_HPP
