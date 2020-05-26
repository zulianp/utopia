#ifndef UTOPIA_EVAL_DETERMINANT_HPP
#define UTOPIA_EVAL_DETERMINANT_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    //[implementors guidelines]
    //[minimal] this provides a cross-backen implementation of the determinant computation for n <= 3
    //[optimized] provide a backend specific implementation by specializing the templa

    /// backend generic det for small matrices
    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class EvalDeterminant {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;

        inline static Scalar det_3(const Scalar m00,
                                   const Scalar m01,
                                   const Scalar m02,
                                   const Scalar m10,
                                   const Scalar m11,
                                   const Scalar m12,
                                   const Scalar m20,
                                   const Scalar m21,
                                   const Scalar m22) {
            return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 - m00 * m12 * m21 - m01 * m10 * m22 -
                   m02 * m11 * m20;
        }

        inline static Scalar apply(const Matrix &t) {
            auto s = size(t);
            assert(s.get(0) == s.get(1));

            Read<Matrix> r(t);
            Scalar out;

            switch (s.get(0)) {
                case 1: {
                    out = t.get(0, 0);
                    break;
                }
                case 2: {
                    out = (t.get(0, 0) * t.get(1, 1) - t.get(1, 0) * t.get(0, 1));
                    break;
                }
                case 3: {
                    out = det_3(t.get(0, 0),
                                t.get(0, 1),
                                t.get(0, 2),
                                t.get(1, 0),
                                t.get(1, 1),
                                t.get(1, 2),
                                t.get(2, 0),
                                t.get(2, 1),
                                t.get(2, 2));
                    break;
                }
                case 4: {
                    const Scalar m00 = t.get(0, 0);
                    const Scalar m01 = t.get(0, 1);
                    const Scalar m02 = t.get(0, 2);
                    const Scalar m03 = t.get(0, 3);

                    const Scalar m10 = t.get(1, 0);
                    const Scalar m11 = t.get(1, 1);
                    const Scalar m12 = t.get(1, 2);
                    const Scalar m13 = t.get(1, 3);

                    const Scalar m20 = t.get(2, 0);
                    const Scalar m21 = t.get(2, 1);
                    const Scalar m22 = t.get(2, 2);
                    const Scalar m23 = t.get(2, 3);

                    const Scalar m30 = t.get(3, 0);
                    const Scalar m31 = t.get(3, 1);
                    const Scalar m32 = t.get(3, 2);
                    const Scalar m33 = t.get(3, 3);

                    // https://www.math10.com/en/algebra/matrices/determinant.html
                    out = m00 * det_3(m11, m12, m13, m21, m22, m23, m31, m32, m33)

                          - m01 * det_3(m10, m12, m13, m20, m22, m23, m30, m32, m33)

                          + m02 * det_3(m10, m11, m13, m20, m21, m23, m30, m31, m33)

                          - m03 * det_3(m10, m11, m12, m20, m21, m22, m30, m31, m32);
                    break;
                }
                default: {
                    assert(false && "not implemented");
                    std::cerr << "det not implemented for matrices with n > 3" << std::endl;
                    return -1;
                }
            }

            return out;
        }
    };

    template <class T, class Traits, int Backend>
    class Eval<Determinant<Tensor<T, 2> >, Traits, Backend> {
    public:
        using T2 = utopia::Tensor<T, 2>;
        using Scalar = typename Traits::Scalar;

        template <typename Real>
        inline static void apply(const Determinant<T2> &expr, Number<Real> &num) {
            num = apply(expr);
        }

        inline static Scalar apply(const Determinant<T2> &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            const Scalar ret = EvalDeterminant<T>::apply(expr.expr().derived());
            UTOPIA_TRACE_END(expr);
            return ret;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_EVAL_DETERMINANT_HPP
