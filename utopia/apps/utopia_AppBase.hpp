#ifndef UTOPIA_APP_BASE_HPP
#define UTOPIA_APP_BASE_HPP

#include "utopia_Base.hpp"
#ifdef UTOPIA_ENABLE_VC
#define USE_SIMD_ASSEMBLY
#endif  // UTOPIA_ENABLE_VC

#ifdef USE_SIMD_ASSEMBLY
// #include "utopia_simd_Assembler.hpp"
#include "utopia_simd_Assembler_v2.hpp"
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
