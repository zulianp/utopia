#ifndef UTOPIA_TENSOR_HPP
#define UTOPIA_TENSOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {

	//CRTP type
	template<typename Derived, int Order_>
	class Tensor {
	public:
		static const int Order   = Order_;
		static const int StoreAs = UTOPIA_BY_REFERENCE;

		CONST_DERIVED_CRT(Derived);
		DERIVED_CRT(Derived);
	};

	template<typename Derived, int Order_>
	class Traits<Tensor<Derived, Order_>> : public Traits<Derived> {
	public:
	    static const int Order = Order;
	};

}

#endif //UTOPIA_TENSOR_HPP
