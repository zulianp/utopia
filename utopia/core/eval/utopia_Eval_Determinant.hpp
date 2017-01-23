#ifndef UTOPIA_EVAL_DETERMINANT_HPP
#define UTOPIA_EVAL_DETERMINANT_HPP  

#include "utopia_Eval_Empty.hpp"

namespace utopia {
	///backend generic det for small matrices
	template<class Tensor, class Traits, int Backend>
	class Eval<Determinant< Wrapper<Tensor, 2> >, Traits, Backend> {
	public:
		typedef utopia::Wrapper<Tensor, 2> Tensor2;
		typedef typename Traits::Scalar Scalar;

		inline static Scalar apply(const Determinant<Tensor2> &expr)
		{
			auto &t = expr.expr();
			auto s = size(t);
			assert(s.get(0) == s.get(1));

			Read<Tensor2> r(t);
			switch(s.get(0)) {
				case 1:
				{
					return t.get(0, 0);
				}
				case 2: 
				{
					return ( t.get(0, 0) * t.get(1, 1) - t.get(1, 0) * t.get(0, 1) );
				}
				case 3: 
				{
					const Scalar m00 = t.get(0, 0);
					const Scalar m01 = t.get(0, 1);
					const Scalar m02 = t.get(0, 2);
					const Scalar m10 = t.get(1, 0);
					const Scalar m11 = t.get(1, 1);
					const Scalar m12 = t.get(1, 2);
					const Scalar m20 = t.get(2, 0);
					const Scalar m21 = t.get(2, 1);
					const Scalar m22 = t.get(2, 2);

					return m00 * m11 * m22  + 
						   m01 * m12 * m20  + 
						   m02 * m10 * m21  - 
						   m00 * m12 * m21  - 
						   m01 * m10 * m22  - 
						   m02 * m11 * m20; 
				} default: {
					assert(false && "not implemented");
					std::cerr << "det not implemented for matrices with n > 3" << std::endl;
					return -1;
				}
			}
		}
	};
}

#endif //UTOPIA_EVAL_DETERMINANT_HPP
