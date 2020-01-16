#ifndef UTOPIA_PETSC_EVAL_KRONECKERPRODUCT_HPP
#define UTOPIA_PETSC_EVAL_KRONECKERPRODUCT_HPP


#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Eval_KroneckerProduct.hpp"

namespace utopia {

	template<class Matrix, class Vector> 
	class EvalKroneckerProduct<Matrix, Vector, PETSC> {
	public:
		static void apply(const Vector &left, const Vector &right, Matrix &result);
	};
}

#endif //UTOPIA_PETSC_EVAL_KRONECKERPRODUCT_HPP
