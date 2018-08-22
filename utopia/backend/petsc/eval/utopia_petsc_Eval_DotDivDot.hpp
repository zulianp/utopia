#ifndef UTOPIA_PETSC_EVAL_DOT_DIV_DOT_HPP
#define UTOPIA_PETSC_EVAL_DOT_DIV_DOT_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

//Divides
// |	Reduce
// |	|	EMultiplies
// |	|	|	Vector
// |	|	|	Vector
// |	|	Plus
// |	Reduce
// |	|	EMultiplies
// |	|	|	Multiply
// |	|	|	|	SparseMatrix
// |	|	|	|	Vector
// |	|	|	Vector
// |	|	Plus

namespace utopia {
	template<class X, class FunOfX, class Traits>
	class Eval< Binary<
					Dot<Wrapper<X, 1>, Wrapper<X, 1> >,
					Dot<FunOfX, Wrapper<X, 1> >,
					Divides>, Traits, PETSC> {
	public:
		typedef Binary< Dot<Wrapper<X, 1>, Wrapper<X, 1> >,
					    Dot<FunOfX, Wrapper<X, 1> >,
					    Divides> Expr;

	    typedef typename Traits::Scalar Scalar;

	    inline static Scalar apply(const Expr &expr) {
			UTOPIA_TRACE_BEGIN(expr);

			const auto &x1 = expr.left().expr().left().implementation();
			const auto &x2 = expr.left().expr().right().implementation();

			auto &&x3 = Eval<FunOfX, Traits>::apply(expr.right().expr().left());
			const auto &x4 = expr.right().expr().right().implementation();


			PetscScalar result_num = 0., result_denom = 0.;

			PetscErrorCode ierr = 0;
			if(x1.implementation() == x4.implementation()) {
				Vec vecs[2] = { x2.implementation(), x3.implementation() };
				PetscScalar vals[2] = { 0., 0.};
				ierr = VecMDot(x1.implementation(), 2, vecs, vals); assert(ierr == 0);
				result_num = vals[0];
				result_denom = vals[1];
			} else if(x2.implementation() == x4.implementation()) {
				Vec vecs[2] = { x1.implementation(), x3.implementation() };
				PetscScalar vals[2] = { 0., 0.};
				ierr = VecMDot(x2.implementation(), 2, vecs, vals); assert(ierr == 0);
				result_num = vals[0];
				result_denom = vals[1];
			} else {
				ierr = VecDotBegin(x1.implementation(), x2.implementation(), &result_num); assert(ierr == 0);
				ierr = VecDotBegin(x3.implementation(), x4.implementation(), &result_denom); assert(ierr == 0);

				ierr = VecDotEnd(x1.implementation(), x2.implementation(), &result_num); assert(ierr == 0);
				ierr = VecDotEnd(x3.implementation(), x4.implementation(), &result_denom); assert(ierr == 0);
			}

			Scalar r = result_num/result_denom;
			if(result_denom == 0. || ierr != 0) {
				r = 0.;
			}

			UTOPIA_TRACE_END(expr);
			return r;
	    }
	};
}


#endif //UTOPIA_PETSC_EVAL_DOT_DIV_DOT_HPP
