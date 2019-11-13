#ifndef UTOPIA_DEVICE_INVERSE_HPP
#define UTOPIA_DEVICE_INVERSE_HPP

namespace utopia {

    template<class Expr>
    class DeviceDeterminant {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &A)
        {
            Scalar ret = 0.0;

            const SizeType rows = A.rows();
            const SizeType cols = A.cols();

            if(rows > cols) {
                

            } else {

            }

            return ret;
        }

    private:

        template<class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det(const ExprT &t)
        {
            const auto n = t.rows();
        }
        
        template<class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_1(const Expr &t)
        {
            return t(0, 0);
        }

        template<class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_2(const ExprT &t)
        {
            return t(0, 0) * t(1, 1) - t(1, 0) * t(0, 1);
        }

        UTOPIA_INLINE_FUNCTION static Scalar det_2(
            const Scalar &m00,
            const Scalar &m01,
            const Scalar &m02,
            const Scalar &m10
            )
        {
            return m00 * m11 - m10 * m01;
        }

        template<class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_3(const ExprT &t)
        {
            return det_3(
                m(0, 0),
                m(0, 1),
                m(0, 2),
                m(1, 0),
                m(1, 1),
                m(1, 2),
                m(2, 0),
                m(2, 1),
                m(2, 2)
            );      
        }

        UTOPIA_INLINE_FUNCTION static Scalar det_3(
            const Scalar &m00,
            const Scalar &m01,
            const Scalar &m02,
            const Scalar &m10,
            const Scalar &m11,
            const Scalar &m12,
            const Scalar &m20,
            const Scalar &m21,
            const Scalar &m22)
        {
            return m00 * m11 * m22  +
                   m01 * m12 * m20  +
                   m02 * m10 * m21  -
                   m00 * m12 * m21  -
                   m01 * m10 * m22  -
                   m02 * m11 * m20;
        }

        template<class ExprT>
        UTOPIA_INLINE_FUNCTION static Scalar det_4(const ExprT &t)
        {
            return det_4(
                //Row 0
                t(0, 0),
                t(0, 1),
                t(0, 2),
                t(0, 3),
                //Row 1
                t(1, 0),
                t(1, 1),
                t(1, 2),
                t(1, 3),
                //Row 2
                t(2, 0),
                t(2, 1),
                t(2, 2),
                t(2, 3),
                //Row 3
                t(3, 0),
                t(3, 1),
                t(3, 2),
                t(3, 3)
            );

        }

        //https://www.math10.com/en/algebra/matrices/determinant.html
        UTOPIA_INLINE_FUNCTION static Scalar det_4(
            //Row 0
            const Scalar m00,
            const Scalar m01,
            const Scalar m02,
            const Scalar m03,
            //Row 1
            const Scalar m10,
            const Scalar m11,
            const Scalar m12,
            const Scalar m13,
            //Row 2
            const Scalar m20,
            const Scalar m21,
            const Scalar m22,
            const Scalar m23,
            //Row 3
            const Scalar m30,
            const Scalar m31,
            const Scalar m32,
            const Scalar m33
            )
        {
            return
                m00 * det_3(m11, m12, m13,
                            m21, m22, m23,
                            m31, m32, m33)

                - m01 * det_3(m10, m12, m13,
                              m20, m22, m23,
                              m30, m32, m33)

                + m02 * det_3(m10, m11, m13,
                              m20, m21, m23,
                              m30, m31, m33)

                - m03 * det_3(m10, m11, m12,
                              m20, m21, m22,
                              m30, m31, m32);
        }

    };

    template<class Expr>
    class DeviceInverse {
    public:
        using Scalar   = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(Left &left, const Right &right)
        {

            return ret;
        }
    };

}

#endif //UTOPIA_DEVICE_INVERSE_HPP
