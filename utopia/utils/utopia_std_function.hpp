#ifndef UTOPIA_STD_FUNCTION_HPP
#define UTOPIA_STD_FUNCTION_HPP

#include "utopia_Base.hpp"


#ifdef KOKKOS_ENABLE_CUDA

#include <nvfunctional>

namespace utopia {

    template<typename T>
    using function = nvstd::function<T>;

}

#else

#include <functional>

namespace utopia {

    template<typename T>
    using function = std::function<T>;

}

#endif

#endif //UTOPIA_STD_FUNCTION_HPP
