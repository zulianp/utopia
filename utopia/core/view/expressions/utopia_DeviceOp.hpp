#ifndef UTOPIA_DEVICE_OP_HPP
#define UTOPIA_DEVICE_OP_HPP

#include "utopia_Operators.hpp"

namespace utopia {

    template<class Scalar, class Op> class DeviceOp {};

    template<class Scalar>
    class DeviceOp<Scalar, Max> {
    public:
        UTOPIA_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            return (a > b)? a : b;
        }
    };

    template<class Scalar>
    class DeviceOp<Scalar, Min> {
    public:
        UTOPIA_INLINE_FUNCTION Scalar apply(const Scalar &a, const Scalar &b) const {
            return (a < b)? a : b;
        }
    };

    template<class Scalar>
    class DeviceOp<Scalar, Plus> {
    public:
        UTOPIA_INLINE_FUNCTION static Scalar apply(const Scalar &a, const Scalar &b) {
            return a + b;
        }
    };

}

#endif
