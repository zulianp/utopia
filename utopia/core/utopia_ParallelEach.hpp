#ifndef UTOPIA_PARALLEL_EACH_HPP
#define UTOPIA_PARALLEL_EACH_HPP

#include "utopia_Each.hpp"

namespace utopia {
    template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Traits<Tensor>::FILL_TYPE>
    class ParallelEach : public Each<Tensor, Order, FILL_TYPE> {};

    template<class Tensor, class Fun>
    inline void parallel_each_read(const Tensor &v, Fun fun)
    {
        ParallelEach<Tensor>::apply_read(v, fun);
    }

    template<class Tensor, class Fun>
    inline void parallel_each_write(Tensor &v, Fun fun)
    {
        ParallelEach<Tensor>::apply_write(v, fun);
    }
}

#endif //UTOPIA_PARALLEL_EACH_HPP
