#ifndef UTOPIA_XSIMD_HPP
#define UTOPIA_XSIMD_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_XSIMD

#include <xsimd/xsimd.hpp>

namespace utopia {

    namespace host {
        template <typename T, int N>
        using Batch = xsimd::batch<T, N>;

        using xsimd::load_unaligned;
        using xsimd::store_unaligned;

        using xsimd::load_aligned;
        using xsimd::store_aligned;

    }  // namespace host
}  // namespace utopia

#else  // UTOPIA_WITH_XSIMD

namespace utopia {}

#endif  // UTOPIA_WITH_XSIMD

#endif  // UTOPIA_XSIMD_HPP
