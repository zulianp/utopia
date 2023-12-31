#ifndef UTOPIA_SIMD_MATRIX_HPP
#define UTOPIA_SIMD_MATRIX_HPP

#include "Vc/Vc"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

#include "utopia_DeviceIdentity.hpp"

#include "utopia_simd_Ops.hpp"
#include "utopia_simd_Vector.hpp"

namespace utopia {
    namespace simd_v2 {

        // template <typename T, int Dim, typename...>
        // class Vector {};

        template <typename T, int Rows, int Cols, class SIMDType = Vc::Vector<T>>
        class Matrix {
        public:
            static const int Lanes = SIMDType::Size;
            static const int N = Rows * Cols;
            static const int Size = N * Lanes;

            using Ops = utopia::simd_v2::Ops<SIMDType>;

            enum { StoreAs = UTOPIA_BY_REFERENCE };

            // static constexpr size_t data_alignment() { return SIMDType::MemoryAlignment; }
            // alignas(data_alignment())

            T data[Size];

            inline static constexpr int rows() { return Rows; }

            inline static constexpr int cols() { return Cols; }

            inline T &ref(const int component, const int lane) {
                assert(component < N);
                assert(lane < Lanes);
                return data[component * Lanes + lane];
            }

            inline constexpr const T &get(const int component, const int lane) const {
                assert(component < N);
                assert(lane < Lanes);

                return data[component * Lanes + lane];
            }

            inline constexpr const T &get(const int component_i, const int component_j, const int lane) const {
                return (*this)(component_i * Rows + component_j, lane);
            }

            inline constexpr const T *block(const int idx) const {
                assert(idx < N);
                assert(idx >= 0);

                return &data[idx * Lanes];
            }

            inline T *block(const int i, const int j) {
                assert(i < Rows);
                assert(j < Cols);
                assert(i >= 0);
                assert(j >= 0);

                return block(i * Rows + j);
            }

            inline constexpr const T *block(const int i, const int j) const {
                assert(i < Rows);
                assert(j < Cols);
                assert(i >= 0);
                assert(j >= 0);

                return block(i * Rows + j);
            }

            inline T *block(const int idx) {
                assert(idx < N);
                assert(idx >= 0);

                return &data[idx * Lanes];
            }

            inline SIMDType load(const int i, const int j) const { return Ops::construct(block(i, j)); }

            inline CONST_SIMD_RET SIMDType operator()(const int i, const int j) const {
                return Ops::construct(block(i, j));
            }

            void set(const T &val) {
                for (int i = 0; i < Size; ++i) {
                    data[i] = val;
                }
            }

            void set(const int i, const int j, const SIMDType &val) { Ops::store(val, block(i, j)); }

            inline void dot(const Matrix &other, SIMDType &result) const {
                Ops::template dot<N>(data, other.data, result);
            }

            friend inline constexpr SIMDType dot(const Matrix &l, const Matrix &r) {
                SIMDType ret;
                l.dot(r, ret);
                return ret;
            }

            friend inline constexpr SIMDType inner(const Matrix &l, const Matrix &r) {
                SIMDType ret;
                l.dot(r, ret);
                return ret;
            }

            template <typename Derived>
            friend inline constexpr SIMDType inner(const DeviceExpression<Derived> &l, const Matrix &r) {
                auto &&d = l.derived();
                SIMDType ret = T(0);

                for (int i = 0; i < Rows; ++i) {
                    for (int j = 0; j < Cols; ++j) {
                        ret += d(i, j) * r(i, j);
                    }
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

            friend inline DeviceBinary<DeviceNumber<T>, Matrix, Multiplies> operator*(const Matrix &left,
                                                                                      const T &right) {
                return DeviceBinary<DeviceNumber<T>, Matrix, Multiplies>(right, left);
            }

            friend inline DeviceBinary<DeviceNumber<T>, Matrix, Multiplies> operator*(const T &left,
                                                                                      const Matrix &right) {
                return DeviceBinary<DeviceNumber<T>, Matrix, Multiplies>(left, right);
            }

            friend inline DeviceBinary<DeviceNumber<SIMDType>, Matrix, Multiplies> operator*(const SIMDType &left,
                                                                                             const Matrix &right) {
                return DeviceBinary<DeviceNumber<SIMDType>, Matrix, Multiplies>(left, right);
            }

            friend inline DeviceMultiply<Matrix, Vector<T, Cols, SIMDType>> operator*(
                const Matrix &left,
                const Vector<T, Cols, SIMDType> &right) {
                return DeviceMultiply<Matrix, Vector<T, Cols, SIMDType>>(left, right);
            }

            template <class Expr>
            inline Matrix &operator=(const DeviceExpression<Expr> &expr) {
                auto &&d = expr.derived();
                for (int i = 0; i < Rows; ++i) {
                    for (int j = 0; j < Cols; ++j) {
                        Ops::store(d(i, j), block(i, j));
                    }
                }

                return *this;
            }

            template <class Expr>
            inline Matrix &operator+=(const DeviceExpression<Expr> &expr) {
                auto &&d = expr.derived();
                for (int i = 0; i < Rows; ++i) {
                    for (int j = 0; j < Cols; ++j) {
                        Ops::in_place_add(d(i, j), block(i, j));
                    }
                }

                return *this;
            }

            inline Matrix &operator+=(const Matrix &expr) {
                for (int i = 0; i < N; ++i) {
                    auto block_i = block(i);
                    Ops::store(Ops::load(block_i), Ops::load(expr.block(i)), block_i);
                }
                return *this;
            }

            template <typename Factor>
            inline void scale(const Factor &factor) {
                for (int i = 0; i < N; ++i) {
                    Ops::scale(factor, block(i));
                }
            }

            inline Matrix &operator*=(const SIMDType &factor) {
                scale(factor);
                return *this;
            }

            inline Matrix &operator*=(const T &factor) {
                scale(factor);
                return *this;
            }

            friend inline constexpr SIMDType trace(const Matrix &m) {
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
                        Ops::add_scale_store(0.5, block(i, j), block(j, i));
                    }
                }
            }

            inline std::string get_class() const { return "simd::Matrix"; }
        };
    }  // namespace simd_v2

    template <typename T, int Rows, int Cols, typename SIMDType_>
    class Traits<simd_v2::Matrix<T, Rows, Cols, SIMDType_>> {
    public:
        using SizeType = int;
        using Primitive = T;
        using Scalar = SIMDType_;
        static const int Order = 2;
    };

    template <typename T, class Right>
    inline DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies> operator*(
        const Vc::Vector<T> &left,
        const DeviceExpression<Right> &right) {
        return DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies>(left, right.derived());
    }

}  // namespace utopia

#endif  // UTOPIA_SIMD_MATRIX_HPP
