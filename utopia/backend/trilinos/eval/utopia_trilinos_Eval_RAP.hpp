#ifndef UTOPIA_TRILINOS_EVAL_RAP_HPP
#define UTOPIA_TRILINOS_EVAL_RAP_HPP

#include "utopia_trilinos_Types.hpp"

//FIXME find right macro
#ifdef WITH_TPETRAEXT

#include "TpetraExt_TripleMatrixMultiply_decl.hpp"

//useful links:
//https://trilinos.org/docs/dev/packages/tpetra/doc/html/namespaceTpetra_1_1TripleMatrixMultiply.html
//see MultiplyRAP

namespace utopia {

	// mat-mat-mat multiplication
	template<class M1, class M2, class M3, class Traits>
	class Eval< Multiply< Multiply< Transposed<Wrapper<M1, 2>>, Wrapper<M2, 2> >, Wrapper<M3, 2> >, Traits, TRILINOS> {
	public:
		typedef utopia::Multiply< Multiply< Transposed<Wrapper<M1, 2>>, Wrapper<M2, 2> >, Wrapper<M3, 2> > Expr;
		typedef EXPR_TYPE(Traits, Expr) Result;

		inline static Result apply(const Expr &expr)
		{
			Result result;

			UTOPIA_TRACE_BEGIN(expr);

			auto &R = expr.left().left().expr();
			auto &A = expr.left().right();
			auto &P = expr.right();

			result.implementation_ptr().reset(
				new typename Result::crs_mat_type(
					raw_type(R)->getDomainMap(),
					0,
					Tpetra::DynamicProfile
				)
			);

			//Performs optimal triple product
			//Ac = R*A*P,
			Tpetra::TripleMatrixMultiply::MultiplyRAP(
				*raw_type(R),
		        true, //transposeR
		        *raw_type(A),
		        false, //transposeA
		        *raw_type(P),
		        false, //transposeP
		        result.implementation(),
		        true  //call_FillComplete_on_result
			);

			// result.set_domain_and_range(raw_type(P)->getDomainMap(), raw_type(R)->getRangeMap());
			// result.finalize();

			UTOPIA_TRACE_END(expr);
				// assert(result.same_type(Eval<Wrapper<M3, 2>, Traits>::apply(expr.right())));
			return result;
		}
	};

}

#endif //WITH_TPETRAEXT

#endif //UTOPIA_TRILINOS_EVAL_RAP_HPP

