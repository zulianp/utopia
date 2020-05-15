//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_MUTABLE_HPP
#define UTOPIA_UTOPIA_MUTABLE_HPP

#include "utopia_Operators.hpp"

namespace utopia {

    template <class Implementation, class Derived>
    class Mutable {
    public:
        using Scalar = typename Traits<Implementation>::Scalar;

        template <class OtherDerived>
        Derived &operator*=(const Expression<OtherDerived> &expr) {
            derived().evaluator().eval(InPlace<Derived, OtherDerived, Multiplies>(derived(), expr.derived()));
            return derived();
        }

        template <class OtherDerived>
        Derived &operator+=(const Expression<OtherDerived> &expr) {
            derived().evaluator().eval(InPlace<Derived, OtherDerived, Plus>(derived(), expr.derived()));
            return derived();
        }

        template <class OtherDerived>
        Derived &operator-=(const Expression<OtherDerived> &expr) {
            derived().evaluator().eval(InPlace<Derived, OtherDerived, Minus>(derived(), expr.derived()));
            return derived();
        }

        template <class OtherDerived>
        Derived &operator/=(const Expression<OtherDerived> &expr) {
            derived().evaluator().eval(InPlace<Derived, OtherDerived, Divides>(derived(), expr.derived()));
            return derived();
        }

        Derived &operator*=(const Scalar value) {
            derived().evaluator().eval(InPlace<Derived, Number<Scalar>, Multiplies>(derived(), value));
            return derived();
        }

    private:
        inline Derived &derived() { return static_cast<Derived &>(*this); }
        inline const Derived &derived() const { return static_cast<const Derived &>(*this); }
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_MUTABLE_HPP
