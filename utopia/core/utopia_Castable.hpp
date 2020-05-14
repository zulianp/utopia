#ifndef UTOPIA_CASTABLE_HPP
#define UTOPIA_CASTABLE_HPP

#include "utopia_Assign.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template <class Derived, int Order = Derived::Order>
    class Castable {};

    template <class Derived>
    class Castable<Derived, 0> {
    public:
        typedef typename utopia::Traits<Derived>::Scalar Scalar;

        // virtual ~Castable() {}

        inline operator Scalar() const {
            using C = utopia::Construct<Number<Scalar>, Derived>;
            Number<Scalar> ret = 0;
            Eval<C, Traits<Derived>, Traits<Derived>::Backend>::apply(C(ret, derived()));
            return ret;
        }

    private:
        CONST_DERIVED_CRT(Derived)
    };
}  // namespace utopia

#endif  // UTOPIA_CASTABLE_HPP
