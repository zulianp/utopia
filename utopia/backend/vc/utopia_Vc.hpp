#ifndef UTOPIA_VC_HPP
#define UTOPIA_VC_HPP

#include "utopia_Base.hpp"

#include <Vc/Vc>

namespace utopia {

    template <typename T>
    struct ScalarType {
        using Type = T;
    };

    template <typename T>
    struct MemSetZero {
        static void apply(T *data, ptrdiff_t count) { std::memset(data, 0, sizeof(T) * count); }
        inline static constexpr T make() { return static_cast<T>(0); }
    };

    template <typename T>
    struct Disjunction {};

    template <>
    struct Disjunction<bool> {
        UTOPIA_INLINE_FUNCTION constexpr static bool eval(const bool v) { return v; }
    };

    template <typename T>
    struct ScalarType<Vc::Vector<T>> {
        using Type = T;
    };

    template <typename T>
    struct MemSetZero<Vc::Vector<T>> {
        static void apply(Vc::Vector<T> *data, ptrdiff_t count) {
            for (ptrdiff_t i = 0; i < count; ++i) {
                data[i] = Vc::Vector<T>::MemSetZero();
            }
        }

        inline static Vc::Vector<T> make() { return Vc::Vector<T>::MemSetZero(); }
    };

    template <typename... Args>
    struct Disjunction<Vc::Mask<Args...>> {
        inline constexpr static bool eval(const Vc::Mask<Args...> &v) { return v.count() > 0; }
    };

    namespace device {
        template <typename T>
        inline Vc::Vector<T> min(const Vc::Vector<T> &left, const Vc::Vector<T> &right) {
            return std::min(left, right);
        }

        template <typename T>
        inline Vc::Vector<T> max(const Vc::Vector<T> &left, const Vc::Vector<T> &right) {
            return std::max(left, right);
        }

        template <typename T>
        inline Vc::Vector<T> abs(const Vc::Vector<T> &v) {
            return std::abs(v);
        }

        template <typename T>
        UTOPIA_INLINE_FUNCTION constexpr bool disjunction(const T v) {
            return Disjunction<T>::eval(v);
        }

        template <typename T>
        UTOPIA_INLINE_FUNCTION T sum(const T v) {
            return v;
        }

        template <typename T>
        T sum(const Vc::Vector<T> &v) {
            return v.sum();
        }
    }  // namespace device

}  // namespace utopia

#endif  // UTOPIA_VC_HPP
