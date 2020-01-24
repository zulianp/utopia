#ifndef UTOPIA_DEVICE_EIGEN_VECTORS_HPP
#define UTOPIA_DEVICE_EIGEN_VECTORS_HPP

#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    /**
     * Algorithm for 2x2 (https://en.wikipedia.org/wiki/Eigenvalue_algorithm#2%C3%972_matrices)
     * Compute eigen values
     * e2 = (A - lambda_1 I) e1 = (A - lambda_2 I)
     *
     *
     *
     * Algorithm for 3x3 (https://en.wikipedia.org/wiki/Eigenvalue_algorithm#Eigenvectors_of_normal_3%C3%973_matrices)
     * Compute eigen-values Lambda
     * N = diag(1/Lambda) * A for obtaining normal matrix
     * E1 = (A - I)^2 and E2 = ((A-I) * (A+I)) and choose non-zero columns e1 and e2 for two eigen-vectors
     *
     */

    template<class Expr>
    class DeviceEigenVectors :  public DeviceExpression<DeviceEigenVectors<Expr>> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        template<class Vector, class ResultMat>
        UTOPIA_INLINE_FUNCTION static void apply_2(const Expr &mat, const Vector &eigen_values, ResultMat &result)
        {
            result.copy(mat);

            Scalar &u1 = result(0, 0), &u2 = result(1, 0);
            Scalar &v1 = result(0, 1), &v2 = result(1, 1);

            u1 -= eigen_values[1];
            v2 -= eigen_values[0];

            //normalizing eigen vectors

            const Scalar norm_u = device::sqrt(u1 * u1 + u2 * u2);
            const Scalar norm_v = device::sqrt(v1 * v1 + v2 * v2);

            u1 /= norm_u;
            u2 /= norm_u;

            v1 /= norm_v;
            v2 /= norm_v;
        }

        template<class Vector, class ResultMat>
        UTOPIA_INLINE_FUNCTION static void apply_3(const Expr &mat, const Vector &eigen_values, ResultMat &result)
        {
            UTOPIA_DEVICE_ASSERT(!mat.is_alias(result));

            bool is_zero_1 = device::approxeq(eigen_values[0], 0.0, device::epsilon<Scalar>());
            bool is_zero_2 = device::approxeq(eigen_values[1], 0.0, device::epsilon<Scalar>());
            bool is_zero_3 = device::approxeq(eigen_values[2], 0.0, device::epsilon<Scalar>());

            if(is_zero_1 || is_zero_2 || is_zero_3 || mat.is_diagonal(device::epsilon<Scalar>()*100)) {
                result.identity();
                return;
            }

            //expressions (not evaluated)
            auto Am1 = mat - eigen_values[0] * device::identity<Scalar>();
            auto Am2 = mat - eigen_values[1] * device::identity<Scalar>();
            auto Am3 = mat - eigen_values[2] * device::identity<Scalar>();

            const SizeType n = utopia::rows(mat);

            //expressions (not evaluated)
            auto E1 = Am2 * Am3;
            auto E2 = Am3 * Am1;
            auto E3 = Am1 * Am2;

            //lazy evaluation (expensive but no new memory allocs)
            const SizeType j1 = find_non_zero_col(n, E1);
            const SizeType j2 = find_non_zero_col(n, E2);
            const SizeType j3 = find_non_zero_col(n, E3);

            //lazy evaluation (expensive but no new memory allocs)
            copy_col(n, j1, E1, 0, result);
            copy_col(n, j2, E2, 1, result);
            copy_col(n, j3, E3, 2, result);

            normalize_col(n, 0, result);
            normalize_col(n, 1, result);
            normalize_col(n, 2, result);
        }

        template<class Vector, class ResultMat>
        UTOPIA_INLINE_FUNCTION static void apply(const Expr &mat, const Vector &eigen_values, ResultMat &result)
        {
            switch(eigen_values.size()) {
                case 2:
                {
                    apply_2(mat, eigen_values, result);
                    return;
                }
                case 3:
                {
                    apply_3(mat, eigen_values, result);
                    return;
                }
                default:
                {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }
        }

    private:

        template<class MatExpr>
        UTOPIA_INLINE_FUNCTION static SizeType find_non_zero_col(const SizeType &n, const MatExpr &mat)
        {
            // const auto tol = device::epsilon<Scalar>();
            const auto tol = 0.0; //Will this be a problem???

            for(SizeType j = 0; j < n; ++j) {

                bool is_zero = true;
                for(SizeType i = 0; i < n; ++i) {

                    if(!device::approxeq(mat(i, j), 0.0, tol)) {
                        is_zero = false;
                        break;
                    }
                }

                if(!is_zero) {
                    return j;
                }
            }

            UTOPIA_DEVICE_ASSERT(false);
            return 0;
        }

        template<class From, class To>
        UTOPIA_INLINE_FUNCTION static void copy_col(
            const SizeType &n,
            const SizeType &from_col,
            const From &from,
            const SizeType &to_col,
            To &to)
        {
            for(SizeType i = 0; i < n; ++i) {
                to(i, to_col) = from(i, from_col);
            }
        }


        template<class Mat>
        UTOPIA_INLINE_FUNCTION static void normalize_col(const SizeType &n, const SizeType &j, Mat &mat)
        {
            Scalar len = 0;
            for(SizeType i = 0; i < n; ++i) {
                const Scalar v = mat(i, j);
                len += v*v;
            }

            len = device::sqrt(len);

            for(SizeType i = 0; i < n; ++i) {
                mat(i, j) /= len;
            }
        }

    };
}

#endif //UTOPIA_DEVICE_EIGEN_VECTORS_HPP
