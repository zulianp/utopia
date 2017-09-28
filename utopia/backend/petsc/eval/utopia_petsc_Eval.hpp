#ifndef UTOPIA_EVAL_PETSC_HPP
#define UTOPIA_EVAL_PETSC_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Backend.hpp"

/*! @file
* Petsc language extensions
*/

namespace utopia {
	template<class Left, class Right, class Traits>
	class Eval<Construct<Left, LocalDiagBlock<Right> >, Traits, PETSC> {
	public:

		inline static bool apply(const Construct<Left, LocalDiagBlock<Right> > & expr)
		{
			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).build_local_diag_block(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right().expr())
				);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
			return ok;
		}
	};

	template<class Left, class Right, class Traits>
	class Eval<MatrixPtAPProduct<Left, Right>, Traits, PETSC> {
	public:
		inline static EXPR_TYPE(Traits, Left) apply(const MatrixPtAPProduct<Left, Right> &expr)
		{
			EXPR_TYPE(Traits, Left) result;

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).triple_product_PtAP(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right()),
				result
				);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
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

				UTOPIA_LOG_BEGIN(expr);

				//Perform optimal triple product
				const bool ok = UTOPIA_BACKEND(Traits).triple_product(
					Eval<Wrapper<M1, 2 >, Traits>::apply(expr.left().left()),
					Eval<Wrapper<M2, 2 >, Traits>::apply(expr.left().right()),
					Eval<Wrapper<M3, 2 >, Traits>::apply(expr.right()),
					result
					);

				ASSERT(ok);

				UTOPIA_LOG_END(expr);
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

			UTOPIA_LOG_BEGIN(expr);

			//Check if left and right operands are the same object
			if(&expr.left().left().expr() == &expr.right()) {
				//Perform optimal triple product
				const bool ok = UTOPIA_BACKEND(Traits).triple_product_PtAP(
					Eval<M1, Traits>::apply(expr.left().right()),
					Eval<M2, Traits>::apply(expr.right()),
					result
					);

				ASSERT(ok);

			} else {
				//Perform general triple product
				//Maybe map to L^T A R operation if available?
				UTOPIA_BACKEND(Traits).apply(
					Eval<Multiply<Transposed<M1>, M2>, Traits>::apply(expr.left()),
					Eval<M1, Traits>::apply(expr.right()),
					Multiplies(),
					result);
			}

			UTOPIA_LOG_END(expr);
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

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).build_local_redistribute(
				Eval<Left,  Traits>::apply(expr.left()),
				Eval<Right, Traits>::apply(expr.right()),
				result
				);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
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

	        ASSERT(ok);

	        UTOPIA_LOG_END(expr);
	        return result;
	    }
	};


	template<class M, class V1, class V2, class Traits>
	class Eval<
			Binary<Multiply<Wrapper<M, 2>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus>,
			Traits,
			PETSC> {
	public:
		typedef utopia::Binary<Multiply<Wrapper<M, 2>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus> Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).mat_mult_add(
				Eval<Wrapper<M, 2>,  Traits>::apply(expr.left().left()),
				Eval<Wrapper<V1, 1>, Traits>::apply(expr.left().right()),
				Eval<Wrapper<V2, 1>, Traits>::apply(expr.right()),
				result
			);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
			return result;
		}
	};



	template<class M, class V1, class V2, class Traits>
	class Eval<
			Binary<Wrapper<V1, 1>, Multiply<Wrapper<M, 2>, Wrapper<V2, 1>>, Plus>,
			Traits,
			PETSC> {
	public:
		typedef utopia::Binary<Wrapper<V1, 1>, Multiply<Wrapper<M, 2>, Wrapper<V2, 1>>, Plus> Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).mat_mult_add(
				Eval<Wrapper<M, 2>,  Traits>::apply(expr.right().left()),
				Eval<Wrapper<V2, 1>, Traits>::apply(expr.right().right()),
				Eval<Wrapper<V1, 1>, Traits>::apply(expr.left()),
				result
			);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
			return result;
		}
	};



	template<class M, class V1, class V2, class Traits>
	class Eval<
			Binary<Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus>,
			Traits,
			PETSC> {
	public:
		typedef utopia::Binary<Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V1, 1>>, Wrapper<V2, 1>, Plus> Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).mat_multT_add(
				Eval<Wrapper<M, 2>,  Traits>::apply(expr.left().left().expr()),
				Eval<Wrapper<V1, 1>, Traits>::apply(expr.left().right()),
				Eval<Wrapper<V2, 1>, Traits>::apply(expr.right()),
				result
			);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
			return result;
		}
	};



	template<class M, class V1, class V2, class Traits>
	class Eval<
			Binary<Wrapper<V1, 1>, Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V2, 1>>, Plus>,
			Traits,
			PETSC> {
	public:
		typedef utopia::Binary<Wrapper<V1, 1>, Multiply<Transposed<Wrapper<M, 2>>, Wrapper<V2, 1>>, Plus> Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			UTOPIA_LOG_BEGIN(expr);

			const bool ok = UTOPIA_BACKEND(Traits).mat_multT_add(
				Eval<Wrapper<M, 2>,  Traits>::apply(expr.right().left().expr()),
				Eval<Wrapper<V2, 1>, Traits>::apply(expr.right().right()),
				Eval<Wrapper<V1, 1>, Traits>::apply(expr.left()),
				result
			);

			ASSERT(ok);

			UTOPIA_LOG_END(expr);
			return result;
		}
	};

	//for later PetscErrorCode MatGetRowMax(Mat mat,Vec v,PetscInt idx[]) c = min(mat, 1); r = min(mat, 0)
}

#endif //UTOPIA_EVAL_PETSC_HPP
