#ifndef UTOPIA_SCALAR_CAST_HPP
#define UTOPIA_SCALAR_CAST_HPP

#include "utopia_Expressions.hpp"

namespace utopia {

    template <typename T, class Derived>
    T scalar_cast(const Expression<Derived> &expr) {
        // typedef utopia::Traits<Derived> Traits;
        // Evaluator<typename Traits::Vector, Traits::Backend> eval;

        // Number<T> ret = 0;
        // eval.eval( Construct< Number<T>, Derived >( ret, expr.derived() ) );
        // return ret;

        using C = utopia::Construct<Number<T>, Derived>;
        Number<T> ret = 0;
        Eval<C, Traits<Derived>, Traits<Derived>::Backend>::apply(C(ret, expr.derived()));
        return ret;
    }

    template <typename T>
    T scalar_cast(const Factory<Zeros, -1> &expr) {
        assert(expr.size().n_dims() == 0 || expr.size().get(0) == 1);
        return 0.;
    }
}  // namespace utopia

#endif  // UTOPIA_SCALAR_CAST_HPP
