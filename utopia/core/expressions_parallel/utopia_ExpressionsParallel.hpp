#ifndef UTOPIA_EXPRESSIONS_PARALLEL_HPP
#define UTOPIA_EXPRESSIONS_PARALLEL_HPP

/** @defgroup parallel_expressions Parallel
 *  @brief      Expressions specialized for parallel programming.
 */

#include "utopia_LocalDiagBlock.hpp"
#include "utopia_LocalRedistribute.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<typename Tensor, int Order>
    int comm_size(const Wrapper<Tensor, Order> &)
    {
        static_assert(Traits<Wrapper<Tensor, Order>>::Backend < HOMEMADE, "implement me for your backend");
        return 1;
    }

    template<typename Tensor, int Order>
    int comm_rank(const Wrapper<Tensor, Order> &)
    {
        static_assert(Traits<Wrapper<Tensor, Order>>::Backend < HOMEMADE, "implement me for your backend");
        return 0;
    }

    template<typename Tensor, int Order>
    void synchronize(Wrapper<Tensor, Order> &)
    {
        static_assert(Traits<Wrapper<Tensor, Order>>::Backend < HOMEMADE, "implement me for your backend");
    }

}



#endif //UTOPIA_EXPRESSIONS_PARALLEL_HPP
