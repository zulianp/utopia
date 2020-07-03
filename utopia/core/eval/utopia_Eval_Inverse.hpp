#ifndef UTOPIA_EVAL_INVERSE_HPP
#define UTOPIA_EVAL_INVERSE_HPP

#include "utopia_Eval_Empty.hpp"
#include "utopia_Layout.hpp"

namespace utopia {

    template <class Matrix, int Backend = Traits<Matrix>::Backend>
    class EvalInverse {
    public:
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        inline static bool apply(const Matrix &m, Matrix &result) {
            auto s = size(m);
            SizeType n = s.get(0);

            assert(n == SizeType(s.get(1)));
            assert(n <= 3);

            const Scalar d = det(m);

            if (result.rows() != n || result.cols() != n) {
                result.dense(serial_layout(n, n));
            }

            result.set(0.0);

            const bool failed = d == 0.;

            if (failed) {
                std::cerr << "[Error] zero determinant returning zero matrix" << std::endl;
            } else {
                Read<Matrix> r(m);
                Write<Matrix> w(result);
                inverse(m, d, result);
            }

            return !failed;
        }

    private:
        inline static void inverse(const Matrix &m, const Scalar det_m, Matrix &result) {
            auto s = size(m);
            switch (s.get(0)) {
                case 1: {
                    result.set(0, 0, 1.0 / m.get(0, 0));
                    break;
                }
                case 2:
                    inverse_2x2(m, det_m, result);
                    break;
                case 3:
                    inverse_3x3(m, det_m, result);
                    break;
                // case 4: inverse_4x4(m, det_m, result); break;
                default: {
                    assert(false && "not implemented");
                    break;
                }
            }
        }

        static inline void inverse_2x2(const Matrix &m, const Scalar det_m, Matrix &result) {
            result.set(0, 0, m.get(1, 1));
            result.set(1, 1, m.get(0, 0));
            result.set(0, 1, -m.get(0, 1));
            result.set(1, 0, -m.get(1, 0));
            result *= 1.0 / det_m;
        }

        static void inverse_3x3(const Matrix &m, const Scalar det_m, Matrix &result) {
            result.set(0, 0, (m.get(1, 1) * m.get(2, 2) - m.get(1, 2) * m.get(2, 1)) / det_m);
            result.set(0, 1, (m.get(0, 2) * m.get(2, 1) - m.get(0, 1) * m.get(2, 2)) / det_m);
            result.set(0, 2, (m.get(0, 1) * m.get(1, 2) - m.get(0, 2) * m.get(1, 1)) / det_m);
            result.set(1, 0, (m.get(1, 2) * m.get(2, 0) - m.get(1, 0) * m.get(2, 2)) / det_m);
            result.set(1, 1, (m.get(0, 0) * m.get(2, 2) - m.get(0, 2) * m.get(2, 0)) / det_m);
            result.set(1, 2, (m.get(0, 2) * m.get(1, 0) - m.get(0, 0) * m.get(1, 2)) / det_m);
            result.set(2, 0, (m.get(1, 0) * m.get(2, 1) - m.get(1, 1) * m.get(2, 0)) / det_m);
            result.set(2, 1, (m.get(0, 1) * m.get(2, 0) - m.get(0, 0) * m.get(2, 1)) / det_m);
            result.set(2, 2, (m.get(0, 0) * m.get(1, 1) - m.get(0, 1) * m.get(1, 0)) / det_m);
        }

        // template< class ResultExpr >
        // static bool inverse_4x4(const Expr &m, ResultExpr &result)
        // {
        //  Matrix< EntryType, 2, 2, DEFAULT_MATRIX_OPTIONS, MatrixStructure<2, 2> > temp;
        //  Matrix< EntryType, 6, 1, DEFAULT_MATRIX_OPTIONS, MatrixStructure<6, 1> > s, c;

        //  temp(0, 0, m.get(0, 0));
        //  temp(0, 1, m.get(0, 1));
        //  temp(1, 0, m.get(1, 0));
        //  temp(1, 1, m.get(1, 1));
        //  s[0] = temp.determinant();

        //  temp(0, 0, m.get(0, 0));
        //  temp(0, 1, m.get(0, 2));
        //  temp(1, 0, m.get(1, 0));
        //  temp(1, 1, m.get(1, 2));
        //  s[1] = temp.determinant();

        //  temp(0, 0, m.get(0, 0));
        //  temp(0, 1, m.get(0, 3));
        //  temp(1, 0, m.get(1, 0));
        //  temp(1, 1, m.get(1, 3));
        //  s[2] = temp.determinant();

        //  temp(0, 0, m.get(0, 1));
        //  temp(0, 1, m.get(0, 2));
        //  temp(1, 0, m.get(1, 1));
        //  temp(1, 1, m.get(1, 2));
        //  s[3] = temp.determinant();

        //  temp(0, 0, m.get(0, 1));
        //  temp(0, 1, m.get(0, 3));
        //  temp(1, 0, m.get(1, 1));
        //  temp(1, 1, m.get(1, 3));
        //  s[4] = temp.determinant();

        //  temp(0, 0, m.get(0, 2));
        //  temp(0, 1, m.get(0, 3));
        //  temp(1, 0, m.get(1, 2));
        //  temp(1, 1, m.get(1, 3));
        //  s[5] = temp.determinant();

        //  ////////////////////////////////// c

        //  temp(0, 0, m.get(2, 2));
        //  temp(0, 1, m.get(2, 3));
        //  temp(1, 0, m.get(3, 2));
        //  temp(1, 1, m.get(3, 3));
        //  c[5] = temp.determinant();

        //  temp(0, 0, m.get(2, 1));
        //  temp(0, 1, m.get(2, 3));
        //  temp(1, 0, m.get(3, 1));
        //  temp(1, 1, m.get(3, 3));
        //  c[4] = temp.determinant();

        //  temp(0, 0, m.get(2, 1));
        //  temp(0, 1, m.get(2, 2));
        //  temp(1, 0, m.get(3, 1));
        //  temp(1, 1, m.get(3, 2));
        //  c[3] = temp.determinant();

        //  temp(0, 0, m.get(2, 0));
        //  temp(0, 1, m.get(2, 3));
        //  temp(1, 0, m.get(3, 0));
        //  temp(1, 1, m.get(3, 3));
        //  c[2] = temp.determinant();

        //  temp(0, 0, m.get(2, 0));
        //  temp(0, 1, m.get(2, 2));
        //  temp(1, 0, m.get(3, 0));
        //  temp(1, 1, m.get(3, 2));
        //  c[1] = temp.determinant();

        //  temp(0, 0, m.get(2, 0));
        //  temp(0, 1, m.get(2, 1));
        //  temp(1, 0, m.get(3, 0));
        //  temp(1, 1, m.get(3, 1));
        //  c[0] = temp.determinant();

        //  const EntryType det = s[0]*c[5] - s[1]*c[4] + s[2]*c[3] + s[3]*c[2] - s[4]*c[1] + s[5]*c[0];
        //  EXPRESS_ASSERT( det != 0.0 );
        //  if(IsNumericalZero(det))
        //      return false;

        //  // row 0
        //  result.set(0, 0, m.get(1, 1)*c[5] - m.get(1, 2)*c[4] + m.get(1, 3)*c[3]);
        //  result.set(0, 1, -m.get(0, 1)*c[5] + m.get(0, 2)*c[4] - m.get(0, 3)*c[3]);
        //  result.set(0, 2, m.get(3, 1)*s[5] - m.get(3, 2)*s[4] + m.get(3, 3)*s[3];)
        //  result.set(0, 3, -m.get(2, 1)*s[5] + m.get(2, 2)*s[4] - m.get(2, 3)*s[3]);

        //  // row 1
        //  result.set(1, 0, -m.get(1, 0)*c[5] + m.get(1, 2)*c[2] - m.get(1, 3)*c[1]);
        //  result.set(1, 1, m.get(0, 0)*c[5] - m.get(0, 2)*c[2] - m.get(0, 3)*c[1]);
        //  result.set(1, 2, -m.get(3, 0)*s[5] + m.get(3, 2)*s[2] - m.get(3, 3)*s[1];)
        //  result.set(1, 3, m.get(2, 0)*s[5] - m.get(2, 2)*s[2] + m.get(2, 3)*s[1]);

        //  // row 2
        //  result.set(2, 0, m.get(1, 0)*c[4] - m.get(1, 1)*c[2] + m.get(1, 3)*c[0]);
        //  result.set(2, 1, -m.get(0, 0)*c[4] + m.get(0, 1)*c[2] - m.get(0, 3)*c[0]);
        //  result.set(2, 2, m.get(3, 0)*s[4] - m.get(3, 1)*s[2] + m.get(3, 3)*s[0]);
        //  result.set(2, 3, -m.get(2, 0)*s[4] + m.get(2, 1)*s[2] - m.get(2, 3)*s[0];)

        //  // row 3
        //  result.set(3, 0, -m.get(1, 0)*c[3] + m.get(1, 1)*c[1] - m.get(1, 2)*c[0]);
        //  result.set(3, 1, m.get(0, 0)*c[3] - m.get(0, 1)*c[1] - m.get(0, 2)*c[0]);
        //  result.set(3, 2, -m.get(3, 0)*s[3] + m.get(3, 1)*s[1] - m.get(3, 2)*s[0];)
        //  result.set(3, 3, m.get(2, 0)*s[3] - m.get(2, 1)*s[1] + m.get(2, 2)*s[0]);

        //  result *= 1.0/det;
        //  return true;
        // }
    };

    template <class Derived, class Traits, int Backend>
    class Eval<Inverse<Tensor<Derived, 2> >, Traits, Backend> {
    public:
        using Tensor2 = Derived;
        typedef utopia::Inverse<Tensor<Derived, 2> > Expr;

        using Scalar = typename Traits::Scalar;

        inline static Derived apply(const Expr &expr) {
            Tensor2 ret;
            apply(expr, ret);
            return ret;
        }

        inline static bool apply(const Expr &expr, Tensor2 &result) {
            UTOPIA_TRACE_BEGIN(expr);
            bool ok = EvalInverse<Tensor2>::apply(expr.expr().derived(), result);
            UTOPIA_TRACE_END(expr);
            return ok;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_EVAL_INVERSE_HPP