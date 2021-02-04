#ifndef UTOPIA_SIMD_COMMON_HPP
#define UTOPIA_SIMD_COMMON_HPP

#include "Vc/Vc"

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
    }  // namespace simd

}  // namespace utopia

#endif  // UTOPIA_SIMD_COMMON_HPP
