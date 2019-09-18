#ifndef UTOPIA_KOKKOS_OPERATIONS_HPP
#define UTOPIA_KOKKOS_OPERATIONS_HPP

#include "utopia_Operators.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <iostream>

namespace utopia {

    template<class Scalar, class Op> class KokkosOp {};

    template<class Scalar>
    class KokkosOp<Scalar, Max> {
    public:
        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            return (a > b)? a : b;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Min> {
    public:
        KOKKOS_INLINE_FUNCTION Scalar apply(const Scalar &a, const Scalar &b) const {
            return (a < b)? a : b;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Plus> {
    public:
        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            return a + b;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Divides> {
    public:
        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            return a / b;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, IsNaNOrInf> {
    public:
        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &val) {
            return
                Kokkos::Details::ArithTraits<Scalar>::isNan(val) ||
                Kokkos::Details::ArithTraits<Scalar>::isInf(val);
        }

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            auto apply_a = apply(a);
            auto apply_b = apply(b);
            return (apply_a > apply_b)? apply_a : apply_b;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Exp> {
    public:
        KokkosOp(const Exp &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::exp(value);
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Sqrt> {
    public:
        KokkosOp(const Sqrt &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::sqrt(value);
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Pow> {
    public:
        KokkosOp(const Pow &in) : a_(in.a_) {}

        KOKKOS_INLINE_FUNCTION Scalar apply(const Scalar &value) const {
            return Kokkos::Details::ArithTraits<Scalar>::pow(value, a_);
        }

        Scalar a_;
    };

    template<class Scalar>
    class KokkosOp<Scalar, Reciprocal<Scalar>> {
    public:
        KokkosOp(const Reciprocal<Scalar> &in) : num_(in.numerator()) {}

        KOKKOS_INLINE_FUNCTION Scalar apply(const Scalar &value) const {
            return num_/value;
        }

        Scalar num_;
    };

    template<class Scalar>
    class KokkosOp<Scalar, Pow2> {
    public:
        KokkosOp(const Pow2 &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return value * value;
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Abs> {
    public:
        KokkosOp(const Abs &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::abs(value);
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Log> {
    public:
        KokkosOp(const Log &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::log(value);
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Sin> {
    public:
        KokkosOp(const Sin &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::sin(value);
        }
    };

    template<class Scalar>
    class KokkosOp<Scalar, Cos> {
    public:
        KokkosOp(const Cos &) {}

        KOKKOS_INLINE_FUNCTION static Scalar apply(const Scalar &value) {
            return Kokkos::Details::ArithTraits<Scalar>::cos(value);
        }
    };
}

#endif
