#ifndef UTOPIA_SCALAR_CAST_HPP
#define UTOPIA_SCALAR_CAST_HPP 

#include "utopia_Expressions.hpp"

namespace utopia {
	
	template<typename T, class Derived>
    T scalar_cast(const Expression<Derived> &expr)
    {   
    	typedef utopia::Traits<Derived> Traits;
        Evaluator<typename Traits::Vector, Traits::Backend> eval;

        Number<T> ret = 0;
        eval.eval( Construct< Number<T>, Derived >( ret, expr.derived() ) );
        return ret;
    }
}

#endif //UTOPIA_SCALAR_CAST_HPP
