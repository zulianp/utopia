#ifndef UTOPIA_EVAL_KRONECKER_PRODUCT_HPP
#define UTOPIA_EVAL_KRONECKER_PRODUCT_HPP

#include "utopia_Base.hpp"
#include "utopia_Eval_Empty.hpp"

namespace utopia {

	template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend> 
	class EvalKroneckerProduct {
	public:

		static void apply(const Vector &left, const Vector &right, Matrix &result)
		{
			static_assert(Backend < HOMEMADE, "implement in backend");
		}

	};

}

#endif //UTOPIA_EVAL_KRONECKER_PRODUCT_HPP
