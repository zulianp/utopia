#ifndef UTOPIA_COMPARABLE_HPP
#define UTOPIA_COMPARABLE_HPP

#include "utopia_Traits.hpp"

namespace utopia {
	template<class Derived>
	class Comparable {
	public:
		virtual ~Comparable() {}
		
		using Scalar = typename utopia::Traits<Derived>::Scalar;
		virtual bool equals(const Derived &other, Scalar tol = 0.0) const = 0;

	};
}

#endif //UTOPIA_COMPARABLE_HPP
