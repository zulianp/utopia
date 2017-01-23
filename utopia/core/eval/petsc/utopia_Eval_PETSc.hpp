#ifndef UTOPIA_EVAL_PETSC_HPP
#define UTOPIA_EVAL_PETSC_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_PETScTraits.hpp"
#include "utopia_PETScBackend.hpp"

/*! @file
* Petsc language extensions
*/

namespace utopia {
	template<class Left, class Right, class Traits>
	class Eval<Construct<Left, LocalDiagBlock<Right> >, Traits, PETSC> {
	public:

		inline static bool apply(const Construct<Left, LocalDiagBlock<Right> > & expr)
		{
			const bool ok = UTOPIA_BACKEND(Traits).build_local_diag_block(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right().expr())
				);

			assert(ok);
			return ok;
		}
	};

	template<class Left, class Right, class Traits>
	class Eval<MatrixPtAPProduct<Left, Right>, Traits, PETSC> {
	public:
		inline static EXPR_TYPE(Traits, Left) apply(const MatrixPtAPProduct<Left, Right> &expr)
		{
			EXPR_TYPE(Traits, Left) result;
			const bool ok = UTOPIA_BACKEND(Traits).triple_product_PtAP(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right()),
				result
				);

			assert(ok);
			return result;
		}
	};

	// general mat-mat-mat multiplication
	template<class M1, class M2, class M3, class Traits>
	class Eval< Multiply< Multiply< Wrapper<M1, 2 >, Wrapper<M2, 2 >>, Wrapper<M3, 2 > >,  Traits,  PETSC> 
	{
		public:
			typedef utopia::Multiply< Multiply< Wrapper<M1, 2 >, Wrapper<M2, 2 > >, Wrapper<M3, 2 > > Expr;

			inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
			{
				EXPR_TYPE(Traits, Expr) result;
				//Perform optimal triple product
				const bool ok = UTOPIA_BACKEND(Traits).triple_product(
					Eval<Wrapper<M1, 2 >, Traits>::apply(expr.left().left()),
					Eval<Wrapper<M2, 2 >, Traits>::apply(expr.left().right()),
					Eval<Wrapper<M3, 2 >, Traits>::apply(expr.right()),
					result
					);

				assert(ok);
				return result;
			}
	};

	//! [pattern matching and optimizations]

	/*! 
	* @brief Triple product (m1^T * m2 * m1 := transpose(m1) * m2 * m1) optimization for the petsc backend 
	*/
	template<class M1, class M2, class Traits>
	class Eval<
		//The pattern to match
		Multiply< Multiply<Transposed<M1>, M2>, M1>, 
		//Type information
		Traits, 
		//Restriction to the backend with the PETSC tag.
		PETSC> {
	public:
		typedef utopia::Multiply< Multiply<Transposed<M1>, M2>, M1> Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			//Check if left and right operands are the same object
			if(&expr.left().left().expr() == &expr.right()) {
				//Perform optimal triple product
				const bool ok = UTOPIA_BACKEND(Traits).triple_product_PtAP(
					Eval<M1, Traits>::apply(expr.left().right()),
					Eval<M2, Traits>::apply(expr.right()),
					result
					);

				assert(ok);

			} else {
				//Perform general triple product
				//Maybe map to L^T A R operation if available?
				UTOPIA_BACKEND(Traits).apply(
					Eval<Multiply<Transposed<M1>, M2>, Traits>::apply(expr.left()),
					Eval<M1, Traits>::apply(expr.right()),
					Multiplies(),
					result);
			}

			return result;
		}
	};

	//! [pattern matching and optimizations]

	template<class Left, class Right, class Traits>
	class Eval<LocalRedistribute<Left, Right>, Traits, PETSC> {
	public:

		inline static EXPR_TYPE(Traits, Left) apply(const LocalRedistribute<Left, Right> &expr)
		{
			EXPR_TYPE(Traits, Left) result;
			const bool ok = UTOPIA_BACKEND(Traits).build_local_redistribute(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right()),
				result
				);

			assert(ok);
			return result;
		}
	};



	template<class Left, class Right, typename ScalarT, class Traits>
	class Eval<Binary<
	                    Binary<Number<ScalarT>, Left,  Multiplies>, 
	                    Binary<Number<ScalarT>, Right, Multiplies>, 
	                    Plus
	                 >,
	                 Traits,
	                 PETSC> {
	public:
		typedef utopia::Binary<
	                    	Binary<Number<ScalarT>, Left,  Multiplies>, 
	                    	Binary<Number<ScalarT>, Right, Multiplies>, 
	                    	Plus> Expr;

	    inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    {
	        EXPR_TYPE(Traits, Expr) result;

	        UTOPIA_LOG_BEGIN(expr);

	        const bool ok = UTOPIA_BACKEND(Traits).waxpby(
	        		expr.left().left(),
	        		Eval<Left,  Traits>::apply(expr.left().right()),
	        		expr.right().left(),
	        		Eval<Right, Traits>::apply(expr.right().right()),
	                result
	        );

	        assert(ok);

	        UTOPIA_LOG_END(expr);
	        return result;
	    }
	};
}

#endif //UTOPIA_EVAL_PETSC_HPP
