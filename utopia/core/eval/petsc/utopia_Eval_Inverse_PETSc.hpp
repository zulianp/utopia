#ifndef UTOPIA_EVAL_INVERSE_PETSC_HPP
#define UTOPIA_EVAL_INVERSE_PETSC_HPP 

#include "utopia_Eval_Empty.hpp"
#include "utopia_Inverse.hpp"

namespace utopia {
	template<class Left, class Right, class Traits>
	class Eval< Construct<Left, Inverse<Right> >, Traits, PETSC> {
	public:
		typedef utopia::Construct<Left, Inverse<Right> > Expr;
	    typedef typename TypeAndFill<Traits, Left>::Type Result;

	    inline static void apply(const Expr &expr) {
	        auto & left   = Eval<Left, Traits>::apply(expr.left());
	        auto && right = Eval<Right, Traits>::apply(expr.right().expr());
       		
	    	const bool ok = UTOPIA_BACKEND(Traits).inverse(right, left); assert(ok);
	    }
	};

	template<class Left, class Right, class Traits>
	class Eval< Assign<Left, Inverse<Right> >, Traits, PETSC> {
	public:
		typedef utopia::Assign<Left, Inverse<Right> > Expr;
	    typedef typename TypeAndFill<Traits, Left>::Type Result;

	    inline static void apply(const Expr &expr) {
	        auto & left   = Eval<Left, Traits>::apply(expr.left());
	        auto && right = Eval<Right, Traits>::apply(expr.right().expr());
	      		
	    	const bool ok = UTOPIA_BACKEND(Traits).inverse(right, left); assert(ok);
	    }
	};

}

#endif //UTOPIA_EVAL_INVERSE_PETSC_HPP

