#ifndef UTOPIA_DEVICE_INVERSE_HPP
#define UTOPIA_DEVICE_INVERSE_HPP

#include "utopia_DeviceDeterminant.hpp"

namespace utopia {

    template <class Expr>
    class DeviceInverse : public DeviceExpression<DeviceInverse<Expr>> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION DeviceInverse(const Expr &expr) : expr_(expr) {}

        template <class Result>
        UTOPIA_INLINE_FUNCTION bool apply(Result &result) const {
            return apply(expr_, result);
        }

        template <class Result>
        UTOPIA_INLINE_FUNCTION static bool apply(const Expr &expr, Result &result) {
            const SizeType rows = result.rows();
            const SizeType cols = result.cols();

            if (rows == 4) {
                Scalar d;
                return invert4(expr, result, d);
            }

            const Scalar d = DeviceDeterminant<Expr>::apply(expr);
            if (d == 0.0) return false;

            if (rows != cols) {
                // expr has rows > cols
                const bool row_full_rank = rows < cols;

                switch (rows) {
                    case 1: {
                        return pseudo_inverse<1>(row_full_rank, expr, d, result);
                    }
                    case 2: {
                        return pseudo_inverse<2>(row_full_rank, expr, d, result);
                    }
                    case 3: {
                        return pseudo_inverse<3>(row_full_rank, expr, d, result);
                    }
                    default: {
                        return false;
                    }
                }

            } else {
                inverse(expr, d, result);
            }

            return true;
        }

        template <int Dim, class Result>
        UTOPIA_INLINE_FUNCTION static bool pseudo_inverse(bool row_full_rank,
                                                          const Expr &expr,
                                                          const Scalar &det,
                                                          Result &result) {
            if (row_full_rank) {
                auto mat = transpose(expr) * expr;
                StaticMatrix<Scalar, Dim, Dim> temp;
                inverse(mat, det, temp);
                result = temp * transpose(expr);
                return true;
            } else {
                UTOPIA_DEVICE_ASSERT(false);
                return false;
            }
        }

        template <class ExprT, class Result>
        UTOPIA_INLINE_FUNCTION static void inverse(const ExprT &expr, const Scalar &det, Result &result) {
            const SizeType n = result.rows();

            switch (n) {
                case 1: {
                    result(0, 0) = 1. / det;
                    return;
                }

                case 2: {
                    const Scalar m00 = expr(0, 0);
                    const Scalar m01 = expr(0, 1);
                    const Scalar m10 = expr(1, 0);
                    const Scalar m11 = expr(1, 1);

                    result(0, 0) = m11 / det;
                    result(1, 1) = m00 / det;
                    result(0, 1) = -m01 / det;
                    result(1, 0) = -m10 / det;
                    return;
                }

                case 3: {
                    const Scalar m00 = expr(0, 0);
                    const Scalar m01 = expr(0, 1);
                    const Scalar m02 = expr(0, 2);

                    const Scalar m10 = expr(1, 0);
                    const Scalar m11 = expr(1, 1);
                    const Scalar m12 = expr(1, 2);

                    const Scalar m20 = expr(2, 0);
                    const Scalar m21 = expr(2, 1);
                    const Scalar m22 = expr(2, 2);

                    result(0, 0) = (m11 * m22 - m12 * m21) / det;
                    result(0, 1) = (m02 * m21 - m01 * m22) / det;
                    result(0, 2) = (m01 * m12 - m02 * m11) / det;
                    result(1, 0) = (m12 * m20 - m10 * m22) / det;
                    result(1, 1) = (m00 * m22 - m02 * m20) / det;
                    result(1, 2) = (m02 * m10 - m00 * m12) / det;
                    result(2, 0) = (m10 * m21 - m11 * m20) / det;
                    result(2, 1) = (m01 * m20 - m00 * m21) / det;
                    result(2, 2) = (m00 * m11 - m01 * m10) / det;
                    return;
                }

                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    result.set(0.0);
                }
            }
        }

        template <class ExprT, class Result>
        UTOPIA_INLINE_FUNCTION static bool invert4(const ExprT &m, Result &inv, Scalar &det) {
            const auto m00 = m(0, 0);
            const auto m10 = m(1, 0);
            const auto m20 = m(2, 0);
            const auto m30 = m(3, 0);
            const auto m01 = m(0, 1);
            const auto m11 = m(1, 1);
            const auto m21 = m(2, 1);
            const auto m31 = m(3, 1);
            const auto m02 = m(0, 2);
            const auto m12 = m(1, 2);
            const auto m22 = m(2, 2);
            const auto m32 = m(3, 2);
            const auto m03 = m(0, 3);
            const auto m13 = m(1, 3);
            const auto m23 = m(2, 3);
            const auto m33 = m(3, 3);

            inv(0, 0) = m11 * m22 * m33 - m11 * m23 * m32 - m21 * m12 * m33 + m21 * m13 * m32 + m31 * m12 * m23 -
                        m31 * m13 * m22;

            inv(1, 0) = -m10 * m22 * m33 + m10 * m23 * m32 + m20 * m12 * m33 - m20 * m13 * m32 - m30 * m12 * m23 +
                        m30 * m13 * m22;

            inv(2, 0) = m10 * m21 * m33 - m10 * m23 * m31 - m20 * m11 * m33 + m20 * m13 * m31 + m30 * m11 * m23 -
                        m30 * m13 * m21;

            inv(3, 0) = -m10 * m21 * m32 + m10 * m22 * m31 + m20 * m11 * m32 - m20 * m12 * m31 - m30 * m11 * m22 +
                        m30 * m12 * m21;

            inv(0, 1) = -m01 * m22 * m33 + m01 * m23 * m32 + m21 * m02 * m33 - m21 * m03 * m32 - m31 * m02 * m23 +
                        m31 * m03 * m22;

            inv(1, 1) = m00 * m22 * m33 - m00 * m23 * m32 - m20 * m02 * m33 + m20 * m03 * m32 + m30 * m02 * m23 -
                        m30 * m03 * m22;

            inv(2, 1) = -m00 * m21 * m33 + m00 * m23 * m31 + m20 * m01 * m33 - m20 * m03 * m31 - m30 * m01 * m23 +
                        m30 * m03 * m21;

            inv(3, 1) = m00 * m21 * m32 - m00 * m22 * m31 - m20 * m01 * m32 + m20 * m02 * m31 + m30 * m01 * m22 -
                        m30 * m02 * m21;

            inv(0, 2) = m01 * m12 * m33 - m01 * m13 * m32 - m11 * m02 * m33 + m11 * m03 * m32 + m31 * m02 * m13 -
                        m31 * m03 * m12;

            inv(1, 2) = -m00 * m12 * m33 + m00 * m13 * m32 + m10 * m02 * m33 - m10 * m03 * m32 - m30 * m02 * m13 +
                        m30 * m03 * m12;

            inv(2, 2) = m00 * m11 * m33 - m00 * m13 * m31 - m10 * m01 * m33 + m10 * m03 * m31 + m30 * m01 * m13 -
                        m30 * m03 * m11;

            inv(3, 2) = -m00 * m11 * m32 + m00 * m12 * m31 + m10 * m01 * m32 - m10 * m02 * m31 - m30 * m01 * m12 +
                        m30 * m02 * m11;

            inv(0, 3) = -m01 * m12 * m23 + m01 * m13 * m22 + m11 * m02 * m23 - m11 * m03 * m22 - m21 * m02 * m13 +
                        m21 * m03 * m12;

            inv(1, 3) = m00 * m12 * m23 - m00 * m13 * m22 - m10 * m02 * m23 + m10 * m03 * m22 + m20 * m02 * m13 -
                        m20 * m03 * m12;

            inv(2, 3) = -m00 * m11 * m23 + m00 * m13 * m21 + m10 * m01 * m23 - m10 * m03 * m21 - m20 * m01 * m13 +
                        m20 * m03 * m11;

            inv(3, 3) = m00 * m11 * m22 - m00 * m12 * m21 - m10 * m01 * m22 + m10 * m02 * m21 + m20 * m01 * m12 -
                        m20 * m02 * m11;

            det = m00 * inv(0, 0) + m01 * inv(1, 0) + m02 * inv(2, 0) + m03 * inv(3, 0);

            if (det == 0) return false;

            const Scalar inv_det = 1.0 / det;

            inv *= inv_det;
            return true;
        }

        template <class Array, class ArrayResult>
        UTOPIA_INLINE_FUNCTION static bool invert4_array(const Array &m, ArrayResult &inv, Scalar &det) {
            inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] +
                     m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

            inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] -
                     m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

            inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] +
                     m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

            inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] -
                      m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

            inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] -
                     m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

            inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] +
                     m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

            inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] -
                     m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

            inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] +
                      m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

            inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] +
                     m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

            inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] -
                     m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

            inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] +
                      m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

            inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] -
                      m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

            inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] -
                     m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

            inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] +
                     m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

            inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] -
                      m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

            inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] +
                      m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

            det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

            if (det == 0) return false;

            const Scalar inv_det = 1.0 / det;

            for (int i = 0; i < 16; i++) {
                inv[i] *= inv_det;
            }

            return true;
        }

        UTOPIA_INLINE_FUNCTION const Expr &expr() const { return expr_; }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

}  // namespace utopia

#endif  // UTOPIA_DEVICE_INVERSE_HPP
