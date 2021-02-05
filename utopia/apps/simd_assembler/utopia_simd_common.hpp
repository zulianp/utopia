#ifndef UTOPIA_SIMD_COMMON_HPP
#define UTOPIA_SIMD_COMMON_HPP

#include "Vc/Vc"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace simd {

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
            static const int Size = Rows * Cols;
            T data_[Size];

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

            friend inline constexpr T inner(const Matrix &l, const Matrix &r) { return dot(l, r); }
        };
    }  // namespace simd

    template <typename T, int Dim, typename... Args>
    class Traits<simd::Vector<T, Dim, Args...>> {
    public:
        using Scalar = T;
    };

    namespace simd {
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
        class Vector<T, 2> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            T data_[2] = {simd::Zero<T>::value(), simd::Zero<T>::value()};

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
        };

        template <typename T>
        class Vector<T, 3> final /*: public DeviceExpression<Vector<T, 3>>*/ {
        public:
            T data_[3] = {simd::Zero<T>::value(), simd::Zero<T>::value(), simd::Zero<T>::value()};

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
        };

        template <typename T, int Dim>
        inline T inner(const Vector<T, Dim> &l, const Vector<T, Dim> &r) {
            return dot(l, r);
        }

    }  // namespace simd
}  // namespace utopia

#endif  // UTOPIA_SIMD_COMMON_HPP
