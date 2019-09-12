#ifndef UTOPIA_REDUCIBLE_HPP
#define UTOPIA_REDUCIBLE_HPP

#include "utopia_Operators.hpp"

namespace utopia {

	template<typename T>
	class Reducible {
	public:
		virtual ~Reducible() {}

		// virtual void reduce(const Plus &) = 0;
		// virtual void reduce(const AbsPlus &) = 0;
		// virtual void reduce(const Minus &) = 0;
		// virtual void reduce(const Multiplies &) = 0;
		// virtual void reduce(const And &) = 0;
		virtual T reduce(const Min &) const = 0;
		virtual T reduce(const Max &) const = 0;

		virtual T min() const { return reduce(Min()); }
		virtual T max() const { return reduce(Max()); }

	};
}

#endif //UTOPIA_REDUCIBLE_HPP
