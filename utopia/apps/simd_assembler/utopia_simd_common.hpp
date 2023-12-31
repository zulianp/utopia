#ifndef UTOPIA_SIMD_COMMON_HPP
#define UTOPIA_SIMD_COMMON_HPP

#include "Vc/Vc"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

#include "utopia_DeviceIdentity.hpp"

namespace utopia {
    namespace simd_v1 {

        template <typename T>
        struct ScalarType {
            using Type = T;
        };

        template <typename T>
        struct Zero {
            static void apply(T *data, ptrdiff_t count) { std::memset(data, 0, sizeof(T) * count); }
            inline static constexpr T value() { return static_cast<T>(0); }
        };

        template <typename T>
        struct Disjunction {};

        template <>
        struct Disjunction<bool> {
            inline constexpr static bool eval(const bool v) { return v; }
        };

        template <typename T>
        inline constexpr bool disjunction(const T v) {
            return Disjunction<T>::eval(v);
        }

        template <typename T>
        struct ScalarType<Vc::Vector<T>> {
            using Type = T;
        };

        template <typename T>
        struct Zero<Vc::Vector<T>> {
            static void apply(Vc::Vector<T> *data, ptrdiff_t count) {
                for (ptrdiff_t i = 0; i < count; ++i) {
                    data[i] = Vc::Vector<T>::Zero();
                }
            }

            inline static Vc::Vector<T> value() { return Vc::Vector<T>::Zero(); }
        };

        template <typename... Args>
        struct Disjunction<Vc::Mask<Args...>> {
            inline constexpr static bool eval(const Vc::Mask<Args...> &v) { return v.count() > 0; }
        };

        template <typename T>
        T sum(const Vc::Vector<T> &v) {
            return v.sum();
        }

        template <typename T, int Dim, typename...>
        class Vector {};

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

            void set(const int i, const T &val) { data_[i] = val; }

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

            inline std::string get_class() const { return "simd_v1::Matrix"; }
        };
    }  // namespace simd_v1

    template <typename T, int Dim, typename... Args>
    class Traits<simd_v1::Vector<T, Dim, Args...>> {
    public:
        using Scalar = T;
    };

    template <typename T, int Rows, int Cols>
    class Traits<simd_v1::Matrix<T, Rows, Cols>> {
    public:
        using Scalar = T;
        using SizeType = int;
        static const int Order = 2;
    };

    template <typename T>
    class Traits<Vc::Vector<T>> {
    public:
        using Scalar = Vc::Vector<T>;
        using SizeType = int;
    };

    template <typename T, class Right>
    inline DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies> operator*(
        const Vc::Vector<T> &left,
        const DeviceExpression<Right> &right) {
        return DeviceBinary<DeviceNumber<Vc::Vector<T>>, Right, Multiplies>(left, right.derived());
    }

    namespace simd_v1 {
        template <typename T>
        T sum(const T v) {
            return v;
        }

        template <typename T>
        inline T integrate(const T &v) {
            return sum(v);
        }

        template <typename T>
        inline T integrate(const Vc::Vector<T> &v) {
            return v.sum();
        }

        template <typename T>
        inline Vc::Vector<T> inner(const Vc::Vector<T> &l, const Vc::Vector<T> &r) {
            return l * r;
        }

        template <typename T>
        class Vector<T, 2> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            using SIMDType = T;
            enum { StoreAs = UTOPIA_BY_REFERENCE };

            T data_[2] = {simd_v1::Zero<T>::value(), simd_v1::Zero<T>::value()};

            inline T &x() { return data_[0]; }
            inline constexpr const T &x() const { return data_[0]; }

            inline T &y() { return data_[1]; }
            inline constexpr const T &y() const { return data_[1]; }

            inline constexpr Vector operator*(const T &scale) { return {x() * scale, y() * scale}; }

            inline constexpr const T &operator()(const int idx) const { return data_[idx]; }
            inline T &operator()(const int idx) { return data_[idx]; }

            inline constexpr const T &operator[](const int idx) const { return data_[idx]; }
            inline T &operator[](const int idx) { return data_[idx]; }

            inline constexpr Vector operator+=(const Vector &other) {
                x() += other.x();
                y() += other.y();

                return *this;
            }

            friend inline constexpr T dot(const Vector &l, const Vector &r) { return l.x() * r.x() + l.y() * r.y(); }

            friend void disp(const Vector &v, std::ostream &os = std::cout) { os << v.x() << " " << v.y() << "\n"; }

            void set(const T &val) {
                for (int i = 0; i < 2; ++i) {
                    data_[i] = val;
                }
            }

            void set(const int i, const T &val) { data_[i] = val; }
        };

        template <typename T>
        class Vector<T, 3> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            using SIMDType = T;
            enum { StoreAs = UTOPIA_BY_REFERENCE };

            T data_[3] = {simd_v1::Zero<T>::value(), simd_v1::Zero<T>::value(), simd_v1::Zero<T>::value()};

            inline T &x() { return data_[0]; }
            inline constexpr const T &x() const { return data_[0]; }

            inline T &y() { return data_[1]; }
            inline constexpr const T &y() const { return data_[1]; }

            inline T &z() { return data_[2]; }
            inline constexpr const T &z() const { return data_[2]; }

            inline constexpr Vector operator*(const T &scale) { return {x() * scale, y() * scale, z() * scale}; }

            inline constexpr const T &operator()(const int idx) const { return data_[idx]; }
            inline T &operator()(const int idx) { return data_[idx]; }

            inline constexpr const T &operator[](const int idx) const { return data_[idx]; }
            inline T &operator[](const int idx) { return data_[idx]; }

            inline constexpr Vector operator+=(const Vector &other) {
                x() += other.x();
                y() += other.y();
                z() += other.z();
                return *this;
            }

            friend inline constexpr T dot(const Vector &l, const Vector &r) {
                return l.x() * r.x() + l.y() * r.y() + l.z() * r.z();
            }

            friend void disp(const Vector &v, std::ostream &os = std::cout) {
                os << v.x() << " " << v.y() << " " << v.z() << "\n";
            }

            void set(const T &val) {
                for (int i = 0; i < 3; ++i) {
                    data_[i] = val;
                }
            }

            void set(const int i, const T &val) { data_[i] = val; }
        };

        template <typename T, int Dim>
        inline T inner(const Vector<T, Dim> &l, const Vector<T, Dim> &r) {
            return dot(l, r);
        }

        template <typename Array, typename T2, int Cols>
        inline auto operator*(const TensorView<Array, 2> &mat, const Vector<T2, Cols> &x)
            -> Vector<decltype(mat(0, 0) * T2()), Array::Rows> {
            using RetT = decltype(mat(0, 0) * T2());

            Vector<RetT, Array::Rows> ret;

            for (int i = 0; i < static_cast<int>(Array::Rows); ++i) {
                ret[i] = Zero<RetT>::value();

                for (int j = 0; j < Cols; ++j) {
                    ret[i] += mat(i, j) * x[j];
                }
            }

            return ret;
        }

    }  // namespace simd_v1
}  // namespace utopia

// namespace utopia {
// namespace simd = utopia::simd_v1;
// namespace simd = utopia::simd_v2;
// }  // namespace utopia

#endif  // UTOPIA_SIMD_COMMON_HPP
