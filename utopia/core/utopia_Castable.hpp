#ifndef UTOPIA_CASTABLE_HPP
#define UTOPIA_CASTABLE_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Assign.hpp"

namespace utopia {

    template<class Derived, int Order = Derived::Order>
    class Castable {};

    template<class Derived>
    class Castable<Derived, 0> {
    public:
        typedef typename utopia::Traits<Derived>::Scalar Scalar;

        virtual ~Castable() {}

        inline operator Scalar() const
        {
            typedef utopia::Traits<Derived> Traits;
            Evaluator<typename Traits::Vector, Traits::Backend> eval;

            Number<Scalar> ret = 0;
            eval.eval( Construct< Number<Scalar>, Derived >( ret, derived() ) );
            return ret;
        }

    private:
        CONST_DERIVED_CRT(Derived)
    };
}


#endif //UTOPIA_CASTABLE_HPP

