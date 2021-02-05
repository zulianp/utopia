#ifndef UTOPIA_APP_BASE_HPP
#define UTOPIA_APP_BASE_HPP

#include "utopia_Base.hpp"
#ifdef UTOPIA_WITH_VC
#define USE_SIMD_ASSEMBLY
#endif  // UTOPIA_WITH_VC

#ifdef USE_SIMD_ASSEMBLY
#include "utopia_simd_Assembler.hpp"
#endif  // USE_SIMD_ASSEMBLY

namespace utopia {

#ifndef USE_SIMD_ASSEMBLY
    namespace simd {
        template <typename T>
        inline T integrate(const T v) {
            return v;
        }
    }   // namespace simd
#endif  // USE_SIMD_ASSEMBLY

}  // namespace utopia

#endif  // UTOPIA_APP_BASE_HPP
