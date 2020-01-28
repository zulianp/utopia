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

            Scalar e0 = eigen_values[0];
            Scalar e1 = eigen_values[1];
            Scalar e2 = eigen_values[2];

            Scalar scale = device::max(device::max(device::abs(e0), device::abs(e1)), device::abs(e2));

            if(scale < device::epsilon<Scalar>()*100) {
                result.identity();
                return;
            }

            if(mat.is_diagonal(device::epsilon<Scalar>())) {
                result.identity();
                return;
            }

            if(mat.has_one_nz_per_col(device::epsilon<Scalar>())) {
                result.copy(mat);
                normalize_col(3, 0, result);
                normalize_col(3, 1, result);
                normalize_col(3, 2, result);
                return;
            }

            // static const Scalar tol = device::epsilon<Scalar>()*100;

            bool is_zero_0 = (e0 == 0.0);
            bool is_zero_1 = (e1 == 0.0);
            bool is_zero_2 = (e2 == 0.0);

            bool is_degenerate[3] = { is_zero_0, is_zero_1, is_zero_2 };



            //expressions (not evaluated)
            auto Am0 = mat - e0 * device::identity<Scalar>();
            auto Am1 = mat - e1 * device::identity<Scalar>();
            auto Am2 = mat - e2 * device::identity<Scalar>();

            const SizeType n = utopia::rows(mat);

            //expressions (not evaluated)
            auto E0 = Am1 * Am2;
            auto E1 = Am2 * Am0;
            auto E2 = Am0 * Am1;

            //lazy evaluation (expensive but no new memory allocs)
            if(!is_zero_0) {
                const SizeType j0 = find_non_zero_col(n, E0);
                if(j0 < n) {
                    copy_col(n, j0, E0, 0, result);
                    normalize_col(n, 0, result);
                } else {
                    is_degenerate[0] = true;
                }
            }

            if(!is_zero_1) {
                const SizeType j1 = find_non_zero_col(n, E1);
                if(j1 < n) {
                    copy_col(n, j1, E1, 1, result);
                    normalize_col(n, 1, result);
                } else {
                    is_degenerate[1] = true;
                }
            }

            if(!is_zero_2) {
                const SizeType j2 = find_non_zero_col(n, E2);
                if(j2 < n) {
                    copy_col(n, j2, E2, 2, result);
                    normalize_col(n, 2, result);
                } else {
                    is_degenerate[2] = true;
                }
            }

            int n_degenerates = int(is_degenerate[0]) + int(is_degenerate[1]) + int(is_degenerate[2]);

            if(n_degenerates == 2) {
                gram_schmidt(is_degenerate, result);
                return;
            }

            // UTOPIA_DEVICE_ASSERT(!is_degenerate[0] && !is_degenerate[1] && !is_degenerate[2]);

            StaticVector<Scalar, 3> u, v;
            if(is_zero_0) {
                UTOPIA_DEVICE_ASSERT(e1 != 0.0);
                UTOPIA_DEVICE_ASSERT(e2 != 0.0);

                result.col(1, u);
                result.col(2, v);

                result.set_col(0, cross(u, v));
                normalize_col(n, 0, result);
            }

            if(is_zero_1) {
                UTOPIA_DEVICE_ASSERT(e0 != 0.0);
                UTOPIA_DEVICE_ASSERT(e2 != 0.0);

                result.col(0, u);
                result.col(2, v);

                result.set_col(1, cross(u, v));
                normalize_col(n, 1, result);
            }

            if(is_zero_2) {
                UTOPIA_DEVICE_ASSERT(e0 != 0.0);
                UTOPIA_DEVICE_ASSERT(e1 != 0.0);

                result.col(0, u);
                result.col(1, v);

                result.set_col(2, cross(u, v));
                normalize_col(n, 2, result);
            }
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
            SizeType arg_max = 0;
            Scalar   max_val = 0;

            for(SizeType j = 0; j < n; ++j) {
                Scalar val = 0.0;
                for(SizeType i = 0; i < n; ++i) {
                    const Scalar v = mat(i, j);
                    val += v * v;
                }

                if(max_val < val) {
                    max_val = val;
                    arg_max = j;
                }
            }

            if(max_val == 0.0) {
                //out of bound index
                return n+1;
            } else {
                return arg_max;
            }
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

                UTOPIA_DEVICE_ASSERT( to(i, to_col) == to(i, to_col) );
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

            UTOPIA_DEVICE_ASSERT(len > 0.0);

            len = device::sqrt(len);

            for(SizeType i = 0; i < n; ++i) {
                mat(i, j) /= len;
            }
        }

        template<class Mat>
        UTOPIA_INLINE_FUNCTION static void gram_schmidt(bool is_degenerate[3], Mat &result)
        {
            StaticVector<Scalar, 3> u[3];

            int vec_idx = -1;
            for(int i = 0; i < 3; ++i) {
                if(!is_degenerate[i]) {
                    vec_idx = i;
                    break;
                }
            }

            UTOPIA_DEVICE_ASSERT(vec_idx >= 0);

            result.col(vec_idx, u[0]);
            u[1].copy(u[0]);

            u[1](0) += 1.0;
            u[1] -= u[0] * (dot(u[0], u[1])/dot(u[0], u[0]));
            u[1] /= norm2(u[1]);

            u[2] = cross(u[1], u[1]);
            u[2] /= norm2(u[2]);

            int idx = 0;
            for(int i = 0; i < 2; ++i, ++idx) {
                if(i == vec_idx) {
                    ++idx;
                }

                result.set_col(idx, u[i + 1]);
            }
        }

    };
}

#endif //UTOPIA_DEVICE_EIGEN_VECTORS_HPP
