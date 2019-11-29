#ifndef UTOPIA_DEVICE_EIGEN_VALUES_HPP
#define UTOPIA_DEVICE_EIGEN_VALUES_HPP

#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    template<class Expr>
    class DeviceEigenValues : public DeviceExpression<DeviceEigenValues<Expr>> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        template<class Result>
        UTOPIA_INLINE_FUNCTION static void apply(const Expr &expr, Result &result)
        {
            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            UTOPIA_DEVICE_ASSERT(rows == cols);
            UTOPIA_DEVICE_ASSERT(rows == result.size());

            switch(rows) {
                case 2:
                {
                    apply_2(expr, result);
                    return;
                }
                default:
                {
                    UTOPIA_DEVICE_ASSERT(false);
                    result.set(0.0);
                    return;
                }
            }
        }

    private:

        template<class Result>
        UTOPIA_INLINE_FUNCTION static void apply_2(const Expr &m, Result &result)
        {
            const Scalar a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
            const Scalar discr = (a+d) * (a+d) - 4.0 * (a*d - c*b);
            result(0) = (a + d - device::sqrt(discr))/2.;
            result(1) = (a + d + device::sqrt(discr))/2.;
        }
    };

    template<class Expr>
    class DeviceSingularValues : public DeviceExpression<DeviceSingularValues<Expr>> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        template<class Result>
        UTOPIA_INLINE_FUNCTION static void apply(const Expr &expr, Result &result)
        {
            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            UTOPIA_DEVICE_ASSERT(rows == cols);
            UTOPIA_DEVICE_ASSERT(rows == result.size());

            switch(rows) {
                case 2:
                {
                    apply_2(expr, result);
                    return;
                }
                default:
                {
                    UTOPIA_DEVICE_ASSERT(false);
                    result.set(0.0);
                    return;
                }
            }
        }

    private:

        template<class Result>
        UTOPIA_INLINE_FUNCTION static void apply_2(const Expr &m, Result &result)
        {
            const Scalar a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1);
            const Scalar a2 = a*a;
            const Scalar b2 = b*b;
            const Scalar c2 = c*c;
            const Scalar d2 = d*d;

            const Scalar s1 = a2 + b2 + c2 + d2;
            const Scalar l = (a2 + b2 - c2 - d2);
            const Scalar r = (a*c + b*d);
            const Scalar s2 = device::sqrt(l*l + 4 * r * r);
            result(0) = device::sqrt(0.5 * (s1 + s2));
            result(1) = device::sqrt(0.5 * (s1 - s2));
        }
    };
}

#endif //UTOPIA_DEVICE_EIGEN_VALUES_HPP
