#ifndef UTOPIA_DEVICE_OP_HPP
#define UTOPIA_DEVICE_OP_HPP

#include "utopia_Operators.hpp"

namespace utopia {

    template <class Scalar, class Op>
    class DeviceOp {};

    template <class Scalar>
    class DeviceOp<Scalar, Max> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) {
            return (a > b) ? a : b;
        }
    };


    template <class Scalar>
    class DeviceOp<Scalar, AbsMax> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) {
            auto aa = device::abs(a);
            auto ab = device::abs(b);
            return (aa > ab) ? aa : ab;
        }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Min> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) {
            return (a < b) ? a : b;
        }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Plus> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) { return a + b; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Minus> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) { return a - b; }

        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &x) { return -x; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Multiplies> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) { return a * b; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, EMultiplies> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) { return a * b; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Divides> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a, const Scalar &b) { return a / b; }
    };

    // TODO
    // template<class Scalar>
    // class DeviceOp<Scalar, Sqrt> {
    // public:
    //     UTOPIA_INLINE_FUNCTION static Scalar apply(const Scalar &a) {
    //         return device::sqrt(a);
    //     }
    // };

    template <class Scalar>
    class DeviceOp<Scalar, Pow2> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a) { return a * a; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Abs> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a) { return a < 0.0 ? -a : a; }
    };

    template <class Scalar>
    class DeviceOp<Scalar, Conjugate> {
    public:
        UTOPIA_INLINE_FUNCTION static constexpr Scalar apply(const Scalar &a) {
            // FIXME Missing conj, it has to work on the device
            return a;
        }
    };
}  // namespace utopia

#endif
