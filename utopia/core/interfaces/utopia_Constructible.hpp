#ifndef UTOPIA_CONSTRUCTIBLE_HPP
#define UTOPIA_CONSTRUCTIBLE_HPP

#include "utopia_Size.hpp"

namespace utopia {
	template<typename Scalar_, typename SizeType_, int Order_>
	class Constructible {};

	template<typename Scalar_, typename SizeType_>
	class Constructible<Scalar_, SizeType_, 2> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~Constructible() {}
		virtual void identity(const Size &s) = 0;
		virtual void zeros(const Size &s) { values(s, 0.0); }
		virtual void values(const Size &s, const Scalar val) = 0;
	};
}

#endif //UTOPIA_CONSTRUCTIBLE_HPP
