#ifndef UTOPIA_PETSC_EVAL_MULTIPLY_HPP
#define UTOPIA_PETSC_EVAL_MULTIPLY_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {

	class PetscEvalTripleMatrixProduct {
	public:
		        // static void mat_mult_add(PetscVector &result, const PetscMatrix &m, const PetscVector &right, const PetscVector &left);
		        // static void mat_mult_t_add(PetscVector &result, const PetscMatrix &m, const PetscVector &right, const PetscVector &left);
        static void ptap(PetscMatrix &result, const PetscMatrix &, const PetscMatrix &);
        static void rart(PetscMatrix &result, const PetscMatrix &, const PetscMatrix &);
        static void abc(PetscMatrix &result, const PetscMatrix &, const PetscMatrix &, const PetscMatrix &);
	};

	    // mat-mat-mat multiplication
	    template<class M1, class M2, class M3, class Traits>
	    class Eval< Multiply< Multiply< Tensor<M1, 2>, Tensor<M2, 2> >, Tensor<M3, 2> >, Traits, PETSC> {
	        public:
	            typedef utopia::Multiply< Multiply< Tensor<M1, 2>, Tensor<M2, 2> >, Tensor<M3, 2> > Expr;

	            inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	            {
	                EXPR_TYPE(Traits, Expr) result;

	                UTOPIA_TRACE_BEGIN(expr);

	                //Performs optimal triple product
	                PetscEvalTripleMatrixProduct::abc(
	                    result,
	                    Eval<Tensor<M1, 2>, Traits>::apply(expr.left().left()),
	                    Eval<Tensor<M2, 2>, Traits>::apply(expr.left().right()),
	                    Eval<Tensor<M3, 2>, Traits>::apply(expr.right())
	                    );

	                UTOPIA_TRACE_END(expr);
	                // assert(result.same_type(Eval<Tensor<M3, 2>, Traits>::apply(expr.right())));
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

	            UTOPIA_TRACE_BEGIN(expr);

	            //Check if left and right operands are the same object
	            if(&expr.left().left().expr() == &expr.right()) {
	                //Perform optimal triple product
	                PetscEvalTripleMatrixProduct::ptap(
	                    result,
	                    Eval<M1, Traits>::apply(expr.left().right()),
	                    Eval<M2, Traits>::apply(expr.right())
	                    );

	            } else {
	                //Perform general triple product
	                PetscEvalTripleMatrixProduct::abc(
	                    result,
	                    Eval<Transposed<M1>, Traits>::apply(expr.left().left()),
	                    Eval<M2, Traits>::apply(expr.left().right()),
	                    Eval<M1, Traits>::apply(expr.right())
	                    );
	            }

	            // assert( result.same_type(Eval<M1, Traits>::apply(expr.right())) );
	            UTOPIA_TRACE_END(expr);
	            return result;
	        }
	    };



	    /*!
	    * @brief Triple product (m1 * m2 * m1^T := m1 * m2 * transpose(m1)) optimization for the petsc backend
	    */
	    template<class M1, class M2, class Traits>
	    class Eval<
	        //The pattern to match
	        Multiply< Multiply<M1, M2>, Transposed<M1>>,
	        //Type information
	        Traits,
	        //Restriction to the backend with the PETSC tag.
	        PETSC> {
	    public:
	        typedef utopia::Multiply< Multiply<M1, M2>, Transposed<M1>> Expr;

	        inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	        {
	            EXPR_TYPE(Traits, Expr) result;

	            UTOPIA_TRACE_BEGIN(expr);

	            //Check if left and right operands are the same object
	            if(&expr.left().left() == &expr.right().expr()) {
	                //Perform optimal triple product
	                PetscEvalTripleMatrixProduct::rart(
	                    result,
	                    Eval<M1, Traits>::apply(expr.left().right()),
	                    Eval<M2, Traits>::apply(expr.right().expr())
	                );

	            } else {
	                //Perform general triple product
	                PetscEvalTripleMatrixProduct::abc(
	                    result,
	                    Eval<M1, Traits>::apply(expr.left().left()),
	                    Eval<M2, Traits>::apply(expr.left().right()),
	                    Eval<Transposed<M1>, Traits>::apply(expr.right())
	                );
	            }

	            // assert( result.same_type(Eval<M1, Traits>::apply(expr.right())) );
	            UTOPIA_TRACE_END(expr);
	            return result;
	        }
	    };

	    //! [pattern matching and optimizations]


	    //TODO

	    // template<class Left, class Right, typename ScalarT, class Traits>
	    // class Eval<Binary<
	    //                     Binary<Number<ScalarT>, Left,  Multiplies>,
	    //                     Binary<Number<ScalarT>, Right, Multiplies>,
	    //                     Plus
	    //                  >,
	    //                  Traits,
	    //                  PETSC> {
	    // public:
	    //     typedef utopia::Binary<
	    //                         Binary<Number<ScalarT>, Left,  Multiplies>,
	    //                         Binary<Number<ScalarT>, Right, Multiplies>,
	    //                         Plus> Expr;

	    //     inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    //     {
	    //         EXPR_TYPE(Traits, Expr) result = Eval<Right, Traits>::apply(expr.right().right());

	    //         UTOPIA_TRACE_BEGIN(expr);

	    //         UTOPIA_BACKEND(Traits).axpby(
	    //                 result,
	    //                 expr.left().left(),
	    //                 Eval<Left,  Traits>::apply(expr.left().right()),
	    //                 expr.right().left()
	    //         );

	    //         UTOPIA_TRACE_END(expr);
	    //         return result;
	    //     }
	    // };

	    // template<class M, class V1, class V2, class Traits>
	    // class Eval<
	    //         Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
	    //         Traits,
	    //         PETSC> {
	    // public:
	    //     typedef utopia::Binary<Multiply<Tensor<M, 2>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus> Expr;

	    //     inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    //     {
	    //         EXPR_TYPE(Traits, Expr) result;

	    //         UTOPIA_TRACE_BEGIN(expr);

	    //         UTOPIA_BACKEND(Traits).mat_mult_add(
	    //             result,
	    //             Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left()),
	    //             Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right()),
	    //             Eval<Tensor<V2, 1>, Traits>::apply(expr.right())
	    //         );

	    //         UTOPIA_TRACE_END(expr);
	    //         return result;
	    //     }
	    // };

	    // template<class M, class V1, class V2, class Traits>
	    // class Eval<
	    //         Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus>,
	    //         Traits,
	    //         PETSC> {
	    // public:
	    //     typedef utopia::Binary<Tensor<V1, 1>, Multiply<Tensor<M, 2>, Tensor<V2, 1>>, Plus> Expr;

	    //     inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    //     {
	    //         EXPR_TYPE(Traits, Expr) result;

	    //         UTOPIA_TRACE_BEGIN(expr);

	    //         UTOPIA_BACKEND(Traits).mat_mult_add(
	    //             result,
	    //             Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left()),
	    //             Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
	    //             Eval<Tensor<V1, 1>, Traits>::apply(expr.left())
	    //         );

	    //         UTOPIA_TRACE_END(expr);
	    //         return result;
	    //     }
	    // };

	    // template<class M, class V1, class V2, class Traits>
	    // class Eval<
	    //         Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus>,
	    //         Traits,
	    //         PETSC> {
	    // public:
	    //     typedef utopia::Binary<Multiply<Transposed<Tensor<M, 2>>, Tensor<V1, 1>>, Tensor<V2, 1>, Plus> Expr;

	    //     inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    //     {
	    //         EXPR_TYPE(Traits, Expr) result;

	    //         UTOPIA_TRACE_BEGIN(expr);

	    //         UTOPIA_BACKEND(Traits).mat_mult_t_add(
	    //             result,
	    //             Eval<Tensor<M, 2>,  Traits>::apply(expr.left().left().expr()),
	    //             Eval<Tensor<V1, 1>, Traits>::apply(expr.left().right()),
	    //             Eval<Tensor<V2, 1>, Traits>::apply(expr.right())
	    //         );

	    //         UTOPIA_TRACE_END(expr);
	    //         return result;
	    //     }
	    // };

	    // template<class M, class V1, class V2, class Traits>
	    // class Eval<
	    //         Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus>,
	    //         Traits,
	    //         PETSC> {
	    // public:
	    //     typedef utopia::Binary<Tensor<V1, 1>, Multiply<Transposed<Tensor<M, 2>>, Tensor<V2, 1>>, Plus> Expr;

	    //     inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
	    //     {
	    //         EXPR_TYPE(Traits, Expr) result;

	    //         UTOPIA_TRACE_BEGIN(expr);

	    //         UTOPIA_BACKEND(Traits).mat_mult_t_add(
	    //             result,
	    //             Eval<Tensor<M, 2>,  Traits>::apply(expr.right().left().expr()),
	    //             Eval<Tensor<V2, 1>, Traits>::apply(expr.right().right()),
	    //             Eval<Tensor<V1, 1>, Traits>::apply(expr.left())
	    //         );

	    //         UTOPIA_TRACE_END(expr);
	    //         return result;
	    //     }
	    // };
}

#endif //UTOPIA_PETSC_EVAL_MULTIPLY_HPP
