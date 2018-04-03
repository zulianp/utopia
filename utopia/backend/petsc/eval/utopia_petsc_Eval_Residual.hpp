#ifndef UTOPIA_PETSC_EVAL_RESIDUAL_HPP
#define UTOPIA_PETSC_EVAL_RESIDUAL_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
	template<class A, class X, class B, class Traits>
	class Eval< Binary<Wrapper<B, 1>, Multiply<Wrapper<A, 2>, Wrapper<X, 1>>, Minus>, Traits, PETSC> {
	public:
		typedef utopia::Binary<Wrapper<B, 1>, Multiply<Wrapper<A, 2>, Wrapper<X, 1>>, Minus> Expr;
	    typedef X Result;

	    inline static Result apply(const Expr &expr) {
			UTOPIA_LOG_BEGIN(expr);

			Result r;
			const auto &a = expr.right().left().implementation();
			const auto &x = expr.right().right().implementation();
			const auto &b = expr.left().implementation();

			r.init(x.communicator(), x.type(), x.local_size(), x.size());
			auto ierr = MatResidual(a.implementation(), b.implementation(), x.implementation(), r.implementation()); assert(ierr == 0);
	    	
			UTOPIA_LOG_END(expr);
			return r;
	    }
	};
}

#endif //UTOPIA_PETSC_EVAL_RESIDUAL_HPP
