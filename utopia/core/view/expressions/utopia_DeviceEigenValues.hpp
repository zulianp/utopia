#ifndef UTOPIA_DEVICE_EIGEN_VALUES_HPP
#define UTOPIA_DEVICE_EIGEN_VALUES_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_DeviceIdentity.hpp"

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
                case 3:
                {
                    apply_3(expr, result);
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


        template<class Result>
        UTOPIA_INLINE_FUNCTION static void apply_3(const Expr &m, Result &result)
        {
            //Algorithm from https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
            const Scalar d1  = m(0, 0);
            const Scalar d2  = m(1, 1);
            const Scalar d3  = m(2, 2);

            const Scalar a12 = m(0, 1);
            const Scalar a13 = m(0, 2);
            const Scalar a23 = m(1, 2);

            const Scalar p1 = a12 * a12 + a13 * a13 + a23 * a23;

            // if(device::approxeq(p1, 0., device::epsilon<Scalar>())) {
            if(p1 == 0.) {
                //diagonal matrix
                result(0) = d1;
                result(1) = d2;
                result(2) = d3;
            } else {
                const Scalar q = (d1 + d2 + d3)/3;
                const Scalar d1mq = d1 - q;
                const Scalar d2mq = d2 - q;
                const Scalar d3mq = d3 - q;

                const Scalar p2 = d1mq * d1mq + d2mq * d2mq + d3mq * d3mq + 2.0 * p1;
                const Scalar p = device::sqrt(p2 / 6);

                const Scalar r = det((1.0/p) * (m - q * device::identity<Scalar>()))/2;

                Scalar cos_phi = 0.0;
                Scalar cos_phi_shifted = 0.0;

                if(r <= -1.0) {
                    cos_phi = 0.5; //cos(pi/3)
                    cos_phi_shifted = -1.; //cos(pi/3 + (2*pi/3));
                } else if(r >= 1.0) {
                    cos_phi = 1; //cos(0)
                    cos_phi_shifted = -0.5; //cos(0 + (2*pi/3));
                } else {
                    const Scalar acos_r = device::acos(r)/3.0;
                    cos_phi = device::cos(acos_r);
                    cos_phi_shifted = device::cos(acos_r + 2.0/3.0 * device::pi<Scalar>());
                }

                result(0) = q + 2 * p * cos_phi;
                result(1) = q + 2 * p * cos_phi_shifted;
                result(2) = 3 * q - result(0) - result(1);
            }

            //sorting eigen-values

            if(result(0) > result(2)) {
                device::swap( result(0), result(2) );
            }

            if(result(0) > result(1)) {
                device::swap( result(0), result(1) );
            }

            if(result(1) > result(2)) {
                device::swap( result(1), result(2) );
            }

            // result.sort();
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

            if(result(0) > result(1)) {
                device::swap( result(0), result(1) );
            }

        }

    };
}

#endif //UTOPIA_DEVICE_EIGEN_VALUES_HPP
