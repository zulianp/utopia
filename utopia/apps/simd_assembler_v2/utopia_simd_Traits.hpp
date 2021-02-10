#ifndef UTOPIA_SIMD_TRAITS_HPP
#define UTOPIA_SIMD_TRAITS_HPP

#include "Vc/Vc"
#include "utopia_DeviceExpression.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    namespace simd_v2 {

        template <typename T>
        class Traits<Vc::Vector<T>> {
        public:
            using Scalar = Vc::Vector<T>;
            using SizeType = int;
        };

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
        T sum(const T v) {
            return v;
        }

        template <typename T>
        T sum(const Vc::Vector<T> &v) {
            return v.sum();
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

    }  // namespace simd_v2

}  // namespace utopia

#endif  // UTOPIA_SIMD_TRAITS_HPP
