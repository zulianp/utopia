#ifndef UTOPIA_EPSILON_HPP
#define UTOPIA_EPSILON_HPP

#include "utopia_Base.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#include <Kokkos_ArithTraits.hpp>
#endif  // UTOPIA_WITH_TRILINOS

#include <cmath>

namespace utopia {
    namespace device {

        template <typename T>
        class Epsilon {
        public:
#ifdef KOKKOS_INLINE_FUNCTION
            static inline constexpr T value() { return Kokkos::Details::ArithTraits<T>::epsilon(); }
#else
            static inline T value() { return std::numeric_limits<T>::epsilon(); }
#endif  // KOKKOS_INLINE_FUNCTION
        };

        template <typename T>
        UTOPIA_INLINE_FUNCTION constexpr T epsilon() {
            return Epsilon<T>::value();
        }
    }  // namespace device
}  // namespace utopia

#endif  // UTOPIA_EPSILON_HPP
