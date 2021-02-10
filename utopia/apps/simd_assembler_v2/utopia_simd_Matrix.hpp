#ifndef UTOPIA_SIMD_MATRIX_HPP
#define UTOPIA_SIMD_MATRIX_HPP

#include "Vc/Vc"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

#include "utopia_DeviceIdentity.hpp"

namespace utopia {
    namespace simd_v2 {

        // template <typename T, int Dim, typename...>
        // class Vector {};

        template <typename T, int Rows, int Cols>
        class Matrix {
        public:
            enum { StoreAs = UTOPIA_BY_REFERENCE };

            static const int Size = Rows * Cols;
            T data_[Size];

            inline static constexpr int rows() { return Rows; }

            inline static constexpr int cols() { return Cols; }

            void set(const T &val) {
                for (int i = 0; i < Size; ++i) {
                    data_[i] = val;
                }
            }

            inline T &operator()(const int i, const int j) { return data_[i * Rows + j]; }
            constexpr const T &operator()(const int i, const int j) const { return data_[i * Rows + j]; }

            friend inline constexpr T dot(const Matrix &l, const Matrix &r) {
                T ret = 0.0;
                for (int i = 0; i < Size; ++i) {
                    ret += l.data_[i] * r.data_[i];
                }

                return ret;
            }

            friend inline DeviceTranspose<Matrix> transpose(const Matrix &mat) { return mat; }

            template <typename Derived>
            friend inline DeviceBinary<Derived, Matrix, Plus> operator+(const DeviceExpression<Derived> &l,
                                                                        const Matrix &r) {
                return DeviceBinary<Derived, Matrix, Plus>(l.derived(), r);
            }

            template <typename Derived>
            friend inline DeviceBinary<Derived, Matrix, Plus> operator+(const Matrix &r,
                                                                        const DeviceExpression<Derived> &l) {
                return DeviceBinary<Derived, Matrix, Plus>(l.derived(), r);
            }

            // friend inline DeviceBinary<Number<T>, Matrix, Multiplies> operator*(const Matrix &left, const T &right) {
            //     return DeviceBinary<Number<T>, Matrix, Multiplies>(right, left);
            // }

            // friend inline DeviceBinary<Number<T>, Matrix, Multiplies> operator*(const T &left, const Matrix &right) {
            //     return DeviceBinary<Number<T>, Matrix, Multiplies>(left, right);
            // }

            template <class Expr>
            UTOPIA_INLINE_FUNCTION Matrix &operator=(const DeviceExpression<Expr> &expr) {
                DeviceAssign<Matrix, Expr>::apply(*this, expr.derived());
                return *this;
            }

            UTOPIA_INLINE_FUNCTION Matrix &operator+=(const Matrix &expr) {
                for (int i = 0; i < Size; ++i) {
                    data_[i] += expr.data_[i];
                }
                // for (int i = 0; i < Rows; ++i) {
                //     for (int j = 0; j < Cols; ++j) {
                //         (*this)(i, j) += expr(i, j);
                //     }
                // }
                return *this;
            }

            UTOPIA_INLINE_FUNCTION Matrix &operator*=(const T &factor) {
                for (auto &d : data_) {
                    d *= factor;
                }

                return *this;
            }

            friend inline Matrix operator*(const Matrix &left, const T &right) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline Matrix operator*(const T &right, const Matrix &left) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline constexpr T inner(const Matrix &l, const Matrix &r) { return dot(l, r); }
            friend inline constexpr T trace(const Matrix &m) {
                auto n = std::min(Rows, Cols);

                T ret = 0.0;
                for (int i = 0; i < n; ++i) {
                    ret += m(i, i);
                }

                return ret;
            }

            void symmetrize() {
                assert(Rows == Cols);
                for (int i = 0; i < Rows; ++i) {
                    for (int j = i + 1; j < Cols; ++j) {
                        auto &a_ij = (*this)(i, j);
                        auto &a_ji = (*this)(j, i);
                        const auto val = 0.5 * (a_ij + a_ji);

                        a_ij = val;
                        a_ji = val;
                    }
                }
            }

            inline std::string get_class() const { return "simd::Matrix"; }
        };
    }  // namespace simd_v2

    template <typename T, int Rows, int Cols>
    class Traits<simd_v2::Matrix<T, Rows, Cols>> {
    public:
        using Scalar = T;
        using SizeType = int;
        static const int Order = 2;
    };

    // template <typename T, class Right>
    // inline DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies> operator*(
    //     const Vc::Vector<T> &left,
    //     const DeviceExpression<Right> &right) {
    //     return DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies>(left, right.derived());
    // }

}  // namespace utopia

#endif  // UTOPIA_SIMD_MATRIX_HPP
