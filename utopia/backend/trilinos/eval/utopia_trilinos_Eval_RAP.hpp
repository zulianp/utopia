#ifndef UTOPIA_TRILINOS_EVAL_RAP_HPP
#define UTOPIA_TRILINOS_EVAL_RAP_HPP

#include "utopia_trilinos_Types.hpp"

//FIXME find right macro
#ifdef HAVE_TPETRAEXT

#include "TpetraExt_TripleMatrixMultiply_decl.hpp"

//useful links:
//https://trilinos.org/docs/dev/packages/tpetra/doc/html/namespaceTpetra_1_1TripleMatrixMultiply.html
//see MultiplyRAP

namespace utopia {

	// mat-mat-mat multiplication
	template<class M1, class M2, class M3, class Traits>
	class Eval< Multiply< Multiply< Wrapper<M1, 2>, Wrapper<M2, 2> >, Wrapper<M3, 2> >, Traits, TRILINOS> {
	public:
		typedef utopia::Multiply< Multiply< Wrapper<M1, 2>, Wrapper<M2, 2> >, Wrapper<M3, 2> > Expr;

		inline static EXPR_TYPE(Traits, Expr) apply(const Expr &expr)
		{
			EXPR_TYPE(Traits, Expr) result;

			UTOPIA_TRACE_BEGIN(expr);

			auto &R = Eval<Wrapper<M1, 2>, Traits>::apply(expr.left().left()).implementation();
			auto &A = Eval<Wrapper<M2, 2>, Traits>::apply(expr.left().right()).implementation();
			auto &P = Eval<Wrapper<M3, 2>, Traits>::apply(expr.right()).implementation();

			result.implementation_ptr().reset(
				new decltype(result)(
					R.implementation().getRowMap(),
					0,
					Tpetra::DynamicProfile
				)
			);

			result.set_domain_and_range(P.implementation().getDomainMap(), R.implementation().getRangeMap());

			//Performs optimal triple product
			//Ac = R*A*P,
			Tpetra::TripleMatrixMultiply::MultiplyRAP(
				R,
		        false, //transposeR
		        A,
		        false, //transposeA
		        P,
		        false, //transposeP
		        result.implementation(),
		        false  //call_FillComplete_on_result
			);

			result.finalize();

			UTOPIA_TRACE_END(expr);
				// assert(result.same_type(Eval<Wrapper<M3, 2>, Traits>::apply(expr.right())));
			return result;
		}
	};

}

#endif //HAVE_TPETRAEXT

#endif //UTOPIA_TRILINOS_EVAL_RAP_HPP

