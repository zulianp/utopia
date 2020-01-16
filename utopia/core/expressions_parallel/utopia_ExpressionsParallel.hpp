#ifndef UTOPIA_EXPRESSIONS_PARALLEL_HPP
#define UTOPIA_EXPRESSIONS_PARALLEL_HPP

/** @defgroup parallel_expressions Parallel
 *  @brief      Expressions specialized for parallel programming.
 */

#include "utopia_LocalDiagBlock.hpp"
#include "utopia_LocalRedistribute.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<typename T, int Order>
    int comm_size(const Tensor<T, Order> &)
    {
        static_assert(Traits<Tensor<T, Order>>::Backend < HOMEMADE, "implement me for your backend");
        return 1;
    }

    template<typename T, int Order>
    int comm_rank(const Tensor<T, Order> &)
    {
        static_assert(Traits<Tensor<T, Order>>::Backend < HOMEMADE, "implement me for your backend");
        return 0;
    }

    template<typename T, int Order>
    void synchronize(Tensor<T, Order> &)
    {
        static_assert(Traits<Tensor<T, Order>>::Backend < HOMEMADE, "implement me for your backend");
    }

}



#endif //UTOPIA_EXPRESSIONS_PARALLEL_HPP
