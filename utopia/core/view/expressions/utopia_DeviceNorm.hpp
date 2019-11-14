#ifndef UTOPIA_DEVICE_NORM_HPP
#define UTOPIA_DEVICE_NORM_HPP

#include "utopia_Norm.hpp"

namespace utopia {

    template<typename Scalar, int Order>
    class DeviceNormAux {};

    template<typename Scalar>
    class DeviceNormAux<Scalar, 1> {
    public:
        UTOPIA_INLINE_FUNCTION static void apply(const Scalar &val, Scalar &ret)
        {
            ret += device::abs(val);
        }

        UTOPIA_INLINE_FUNCTION static void finalize(Scalar &) {}
    };

    template<typename Scalar>
    class DeviceNormAux<Scalar, 2> {
    public:
        UTOPIA_INLINE_FUNCTION static void apply(const Scalar &val, Scalar &ret)
        {
            ret += val * val;
        }

        UTOPIA_INLINE_FUNCTION static void finalize(Scalar &val)
        {
            val = device::sqrt(val);
        }
    };

    template<typename Scalar>
    class DeviceNormAux<Scalar, INFINITY_NORM_TAG> {
    public:
        UTOPIA_INLINE_FUNCTION static void apply(const Scalar &val, Scalar &ret)
        {
            ret = device::max(val, ret);
        }

        UTOPIA_INLINE_FUNCTION static void finalize(Scalar &) {}
    };

    template<class Expr, int Order, int TensorOrder = Traits<Expr>::Order>
    class DeviceNorm {};

    template<class Expr, int Order>
    class DeviceNorm<Expr, Order, 1> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr)
        {
            const SizeType n = expr.size();
            Scalar ret = 0.0;

            for(SizeType i = 0; i < n; ++i) {
                DeviceNormAux<Scalar, Order>::apply(expr(i), ret);
            }

            DeviceNormAux<Scalar, Order>::finalize(ret);
            return ret;
        }
    };

    template<class Expr, int Order>
    class DeviceNorm<Expr, Order, 2> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr)
        {
            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            Scalar ret = 0.0;

            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    DeviceNormAux<Scalar, Order>::apply(expr(i, j), ret);
                }
            }

            DeviceNormAux<Scalar, Order>::finalize(ret);
            return ret;
        }
    };

}

#endif //UTOPIA_DEVICE_NORM_HPP
