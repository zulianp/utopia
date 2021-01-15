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
            const Scalar d = DeviceDeterminant<Expr>::apply(expr);
            if (d == 0.0) return false;

            const SizeType rows = result.rows();
            const SizeType cols = result.cols();

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

        UTOPIA_INLINE_FUNCTION const Expr &expr() const { return expr_; }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

}  // namespace utopia

#endif  // UTOPIA_DEVICE_INVERSE_HPP
