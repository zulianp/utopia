#ifndef UTOPIA_DEVICE_DETERMINANT_HPP
#define UTOPIA_DEVICE_DETERMINANT_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_ViewTraits.hpp"

namespace utopia {

    template <class Expr>
    class DeviceDeterminant {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &t) {
            const SizeType rows = t.rows();
            const SizeType cols = t.cols();

            if (rows > cols) {
                return psuedo_det(cols, rows, transpose(t));
            } else if (cols > rows) {
                return psuedo_det(rows, cols, t);
            } else {
                return det(t);
            }
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar psuedo_det(const SizeType &rows, const SizeType &cols, const ExprT &t) {
            switch (rows) {
                case 1: {
                    Scalar ret = 0.0;

                    for (SizeType j = 0; j < cols; ++j) {
                        const Scalar v = t(0, j);
                        ret += v * v;
                    }

                    return device::sqrt(ret);
                }

                case 2: {
                    auto mult = t * transpose(t);

                    const Scalar m01 = mult(0, 1);

                    return device::sqrt(det_2(mult(0, 0), m01, m01, mult(1, 1)));
                }

                case 3: {
                    auto mult = t * transpose(t);

                    const Scalar m01 = mult(0, 1);
                    const Scalar m02 = mult(0, 2);
                    const Scalar m12 = mult(1, 2);

                    return device::sqrt(det_3(mult(0, 0), m01, m02, m01, mult(1, 1), m12, m02, m12, mult(2, 2)));
                }

                case 4: {
                    auto mult = t * transpose(t);

                    const Scalar m01 = mult(0, 1);
                    const Scalar m02 = mult(0, 2);
                    const Scalar m03 = mult(0, 3);

                    const Scalar m12 = mult(1, 2);
                    const Scalar m13 = mult(1, 3);

                    const Scalar m23 = mult(2, 3);

                    return device::sqrt(det_4(mult(0, 0),
                                              m01,
                                              m02,
                                              m03,
                                              m01,
                                              mult(1, 1),
                                              m12,
                                              m13,
                                              m02,
                                              m12,
                                              mult(2, 2),
                                              m23,
                                              m03,
                                              m13,
                                              m23,
                                              mult(3, 3)));
                }

                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det(const ExprT &t) {
            const auto n = t.rows();
            switch (n) {
                case 1: {
                    return det_1(t);
                }
                case 2: {
                    return det_2(t);
                }
                case 3: {
                    return det_3(t);
                }
                case 4: {
                    return det_4(t);
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_1(const ExprT &t) {
            return t(0, 0);
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_2(const ExprT &t) {
            return t(0, 0) * t(1, 1) - t(1, 0) * t(0, 1);
        }

        UTOPIA_INLINE_FUNCTION static Scalar det_2(const Scalar &m00,
                                                   const Scalar &m01,
                                                   const Scalar &m10,
                                                   const Scalar &m11) {
            return m00 * m11 - m10 * m01;
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_3(const ExprT &t) {
            return det_3(t(0, 0), t(0, 1), t(0, 2), t(1, 0), t(1, 1), t(1, 2), t(2, 0), t(2, 1), t(2, 2));
        }

        UTOPIA_INLINE_FUNCTION static Scalar det_3(const Scalar &m00,
                                                   const Scalar &m01,
                                                   const Scalar &m02,
                                                   const Scalar &m10,
                                                   const Scalar &m11,
                                                   const Scalar &m12,
                                                   const Scalar &m20,
                                                   const Scalar &m21,
                                                   const Scalar &m22) {
            return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 - m00 * m12 * m21 - m01 * m10 * m22 -
                   m02 * m11 * m20;
        }

        template <class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_4(const ExprT &t) {
            return det_4(
                // Row 0
                t(0, 0),
                t(0, 1),
                t(0, 2),
                t(0, 3),
                // Row 1
                t(1, 0),
                t(1, 1),
                t(1, 2),
                t(1, 3),
                // Row 2
                t(2, 0),
                t(2, 1),
                t(2, 2),
                t(2, 3),
                // Row 3
                t(3, 0),
                t(3, 1),
                t(3, 2),
                t(3, 3));
        }

        // https://www.math10.com/en/algebra/matrices/determinant.html
        UTOPIA_INLINE_FUNCTION static Scalar det_4(
            // Row 0
            const Scalar m00,
            const Scalar m01,
            const Scalar m02,
            const Scalar m03,
            // Row 1
            const Scalar m10,
            const Scalar m11,
            const Scalar m12,
            const Scalar m13,
            // Row 2
            const Scalar m20,
            const Scalar m21,
            const Scalar m22,
            const Scalar m23,
            // Row 3
            const Scalar m30,
            const Scalar m31,
            const Scalar m32,
            const Scalar m33) {
            return m00 * det_3(m11, m12, m13, m21, m22, m23, m31, m32, m33)

                   - m01 * det_3(m10, m12, m13, m20, m22, m23, m30, m32, m33)

                   + m02 * det_3(m10, m11, m13, m20, m21, m23, m30, m31, m33)

                   - m03 * det_3(m10, m11, m12, m20, m21, m22, m30, m31, m32);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_DEVICE_DETERMINANT_HPP
