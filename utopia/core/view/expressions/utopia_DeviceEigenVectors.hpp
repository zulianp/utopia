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
        using Scalar = typename Traits<Expr>::Scalar;

        template<class Vector, class ResultMat>
        UTOPIA_INLINE_FUNCTION static void apply_2(const Expr &mat, const Vector &eigen_values, ResultMat &result)
        {
            Vector eig_temp;
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
            Vector eig_temp;
            result.copy(mat);



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

    };
}

#endif //UTOPIA_DEVICE_EIGEN_VECTORS_HPP
