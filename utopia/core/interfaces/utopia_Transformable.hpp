#ifndef UTOPIA_TRANSFORMABLE_HPP
#define UTOPIA_TRANSFORMABLE_HPP

#include "utopia_Operators.hpp"

namespace utopia {

	template<typename T>
	class Transformable {
	public:
		virtual ~Transformable() {}

		virtual void transform(const Sqrt &) = 0;
		virtual void transform(const Pow2 &) = 0;
		virtual void transform(const Log &)  = 0;
		virtual void transform(const Exp &)  = 0;
		virtual void transform(const Cos &)  = 0;
		virtual void transform(const Sin &)  = 0;
		virtual void transform(const Abs &)  = 0;
		virtual void transform(const Minus &) = 0;

		virtual void transform(const Pow &p)  = 0;
		virtual void transform(const Reciprocal<T> &f) = 0;

		template<typename Scalar>
		void reciprocal(const Scalar &num) 
		{
			transform(Reciprocal<Scalar>(num));
		}
	};
}

#endif //UTOPIA_TRANSFORMABLE_HPP
