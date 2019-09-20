#ifndef UTOPIA_REDUCIBLE_HPP
#define UTOPIA_REDUCIBLE_HPP

#include "utopia_Operators.hpp"

namespace utopia {

	template<typename T>
	class Reducible {
	public:
		virtual ~Reducible() {}

		virtual T reduce(const Plus &) const = 0;
		virtual T reduce(const Min &)  const = 0;
		virtual T reduce(const Max &)  const = 0;

		virtual T min() const { return reduce(Min()); }
		virtual T max() const { return reduce(Max()); }
		virtual T sum() const { return reduce(Plus()); }

		// virtual T reduce(const AbsPlus &) const = 0;
		// virtual T reduce(const Minus &) const = 0;
		// virtual T reduce(const Multiplies &) const = 0;
		// virtual T reduce(const And &) const = 0;

	};


	template<typename Scalar, typename SizeType>
	class ReducibleMatrix {
	public:
		virtual ~ReducibleMatrix() {}

		virtual SizeType nnz(const Scalar tol = 0.0) const = 0;
		virtual Scalar trace() const  = 0;
	};
}

#endif //UTOPIA_REDUCIBLE_HPP
