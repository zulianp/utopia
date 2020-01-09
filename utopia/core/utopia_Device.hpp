#ifndef UTOPIA_DEVICE_HPP
#define UTOPIA_DEVICE_HPP

#include "utopia_For.hpp"
#include "utopia_Reduction.hpp"

namespace utopia {

    template<int Backend, typename... TParams>
    class Device {
    public:
        template<typename... Args>
        inline static void parallel_for(Args&&... args)
        {
            using ParallelFor = utopia::ParallelFor<Backend>;
            ParallelFor::apply(std::forward<Args>(args)...);
        }

        template<typename... Args>
        inline static void parallel_reduce(Args &&...args)
        {
            using ParallelReduce = utopia::ParallelReduce<Backend>;
            ParallelReduce::apply(std::forward<Args>(args)...);
        }
    };

}

#endif //UTOPIA_DEVICE_HPP
