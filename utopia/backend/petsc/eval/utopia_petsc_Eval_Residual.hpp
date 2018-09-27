#ifndef UTOPIA_PETSC_EVAL_RESIDUAL_HPP
#define UTOPIA_PETSC_EVAL_RESIDUAL_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
	template<class A, class X, class B>
	using PetscMatResidual = Binary<Wrapper<B, 1>, Multiply<Wrapper<A, 2>, Wrapper<X, 1>>, Minus>;

	template<class A, class X, class B, class Traits>
	class Eval<PetscMatResidual<A, X, B>, Traits, PETSC> {
	public:
		typedef utopia::PetscMatResidual<A, X, B> Expr;
	    typedef X Result;

	    inline static Result apply(const Expr &expr) {
			UTOPIA_TRACE_BEGIN(expr);

			Result r;
			const auto &a = expr.right().left().implementation();
			const auto &x = expr.right().right().implementation();
			const auto &b = expr.left().implementation();

			r.init(b.communicator(), b.type(), b.local_size(), b.size());
			PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
			ierr = MatResidual(a.implementation(), b.implementation(), x.implementation(), r.implementation()); assert(ierr == 0);

			UTOPIA_TRACE_END(expr);
			return r;
	    }
	};

	template<class Left, class A, class X, class B, class Traits>
	class Eval<Assign<Left, PetscMatResidual<A, X, B>>, Traits, PETSC> {
	public:
		typedef utopia::PetscMatResidual<A, X, B> Right;
		typedef utopia::Assign<Left, Right> Expr;

	    inline static bool apply(const Expr &assign_expr) {
			UTOPIA_TRACE_BEGIN(assign_expr);
			auto &&expr = assign_expr.right();

			auto &r = assign_expr.left().implementation();
			const auto &a = expr.right().left().implementation();
			const auto &x = expr.right().right().implementation();
			const auto &b = expr.left().implementation();

			if(r.is_null() || r.size() != b.size()) {
				r.repurpose(b.communicator(), b.type(), b.local_size(), b.size());
			}

			PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
			ierr = MatResidual(a.implementation(), b.implementation(), x.implementation(), r.implementation()); assert(ierr == 0);

			UTOPIA_TRACE_END(assign_expr);
			return true;
	    }
	};
}

#endif //UTOPIA_PETSC_EVAL_RESIDUAL_HPP
