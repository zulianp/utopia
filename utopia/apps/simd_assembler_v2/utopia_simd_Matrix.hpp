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

        template <typename T, int Rows, int Cols, int Lanes = Vc::Vector<T>::Size>
        class Matrix {
        public:
            using SIMDType = Vc::Vector<T>;
            static const int N = Rows * Cols;
            static const int Size = N * Lanes;

            enum { StoreAs = UTOPIA_BY_REFERENCE };

            static constexpr size_t data_alignment() { return SIMDType::MemoryAlignment; }
            alignas(data_alignment()) T data[Size];

            inline static constexpr int rows() { return Rows; }

            inline static constexpr int cols() { return Cols; }

            inline T &operator()(const int component, const int lane) {
                assert(component < N);
                assert(lane < Lanes);
                return data[component * Lanes + lane];
            }

            inline constexpr const T &operator()(const int component, const int lane) const {
                assert(component < N);
                assert(lane < Lanes);

                return data[component * Lanes + lane];
            }

            inline constexpr const T &operator()(const int component_i, const int component_j, const int lane) const {
                return (*this)(component_i * Rows + component_j, lane);
            }

            inline constexpr SIMDType load(const int idx) const {
                assert(idx < N);
                assert(idx >= 0);

                return SIMDType(&data[idx * Lanes], Vc::Aligned);
            }

            inline void store(const int idx, const SIMDType &v) {
                assert(idx < N);
                assert(idx >= 0);

                v.store(&data[idx * Lanes], Vc::Aligned);
            }

            inline constexpr SIMDType load(const int i, const int j) const {
                assert(i < Rows);
                assert(j < Cols);
                assert(i >= 0);
                assert(j >= 0);

                return load(i * Rows + j);
            }

            inline void store(const int i, const int j, const SIMDType &v) {
                assert(i < Rows);
                assert(j < Cols);
                assert(i >= 0);
                assert(j >= 0);

                store(i * Rows + j, v);
            }

            void set(const T &val) {
                for (int i = 0; i < Size; ++i) {
                    data[i] = val;
                }
            }

            friend inline constexpr SIMDType dot(const Matrix &l, const Matrix &r) {
                SIMDType ret = T(0);

                for (int i = 0; i < N; ++i) {
                    ret += l.load(i) * r.load(i);
                }

                return ret;
            }

            friend inline constexpr SIMDType inner(const Matrix &l, const Matrix &r) { return dot(l, r); }

            // friend inline DeviceTranspose<Matrix> transpose(const Matrix &mat) { return mat; }

            // template <typename Derived>
            // friend inline DeviceBinary<Derived, Matrix, Plus> operator+(const DeviceExpression<Derived> &l,
            //                                                             const Matrix &r) {
            //     return DeviceBinary<Derived, Matrix, Plus>(l.derived(), r);
            // }

            // template <typename Derived>
            // friend inline DeviceBinary<Derived, Matrix, Plus> operator+(const Matrix &r,
            //                                                             const DeviceExpression<Derived> &l) {
            //     return DeviceBinary<Derived, Matrix, Plus>(l.derived(), r);
            // }

            // friend inline DeviceBinary<Number<T>, Matrix, Multiplies> operator*(const Matrix &left, const T &right) {
            //     return DeviceBinary<Number<T>, Matrix, Multiplies>(right, left);
            // }

            // friend inline DeviceBinary<Number<T>, Matrix, Multiplies> operator*(const T &left, const Matrix &right) {
            //     return DeviceBinary<Number<T>, Matrix, Multiplies>(left, right);
            // }

            // template <class Expr>
            // UTOPIA_INLINE_FUNCTION Matrix &operator=(const DeviceExpression<Expr> &expr) {
            //     DeviceAssign<Matrix, Expr>::apply(*this, expr.derived());
            //     return *this;
            // }

            // UTOPIA_INLINE_FUNCTION Matrix &operator+=(const Matrix &expr) {
            //     for (int i = 0; i < Size; ++i) {
            //         data[i] += expr.data[i];
            //     }
            //     // for (int i = 0; i < Rows; ++i) {
            //     //     for (int j = 0; j < Cols; ++j) {
            //     //         (*this)(i, j) += expr(i, j);
            //     //     }
            //     // }
            //     return *this;
            // }

            template <typename Factor>
            inline void scale(const Factor &factor) {
                for (int i = 0; i < N; ++i) {
                    store(i, load(i) * factor);
                }
            }

            UTOPIA_INLINE_FUNCTION Matrix &operator*=(const T &factor) {
                scale(factor);

                return *this;
            }

            friend inline Matrix operator*(const Matrix &left, const T &right) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline Matrix operator*(const Matrix &left, const SIMDType &right) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline Matrix operator*(const T &right, const Matrix &left) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline Matrix operator*(const SIMDType &right, const Matrix &left) {
                Matrix ret = left;
                ret *= right;
                return ret;
            }

            friend inline constexpr T trace(const Matrix &m) {
                constexpr auto n = std::min(Rows, Cols);

                SIMDType ret = T(0);
                for (int i = 0; i < n; ++i) {
                    ret += m.load(i, i);
                }

                return ret;
            }

            void symmetrize() {
                assert(Rows == Cols);

                for (int i = 0; i < Rows; ++i) {
                    for (int j = i + 1; j < Cols; ++j) {
                        auto a_ij = load(i, j);
                        auto a_ji = load(j, i);

                        const auto val = 0.5 * (a_ij + a_ji);
                        store(i, j, val);
                        store(j, i, val);
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
