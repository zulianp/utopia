#ifndef UTOPIA_DEVICE_HPP
#define UTOPIA_DEVICE_HPP

#include "utopia_For.hpp"

namespace utopia {

    template<int Backend, typename... TParams>
    class Device {
    public:
        template<typename... Args>
        inline static void parallel_for(Args&&... args)
        {
            using ForLoop  = utopia::ParallelFor<Backend>;
            return ForLoop::apply(std::forward<Args>(args)...);
        }
    };

}

#endif //UTOPIA_DEVICE_HPP
