#ifndef UTOPIA_MULTI_REDUCE_HPP
#define UTOPIA_MULTI_REDUCE_HPP

#include "utopia_Traits.hpp"
#include <algorithm>

namespace utopia {

    template<class Tensor, int Backend = Traits<Tensor>::Backend>
    class MultiReduce {
    public:
        using Scalar = typename Traits<Tensor>::Scalar;

        static Scalar multi_min(const Tensor &t1, const Tensor &t2)
        {
            const Scalar min1 = min(t1);
            const Scalar min2 = min(t2);
            return std::min(min1, min2);
        }

    };

    template<class Derived, int Order>
    inline auto multi_min(
        const Tensor<Derived, Order> &t1,
        const Tensor<Derived, Order> &t2) -> typename Traits<Derived>::Scalar
    {
        return MultiReduce<Derived>::multi_min(t1.derived(), t2.derived());
    }
}

#endif //UTOPIA_MULTI_REDUCE_HPP
