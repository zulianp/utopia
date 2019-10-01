#ifndef UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP
#define UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP

#include "utopia_trilinos_ForwardDeclarations.hpp"

namespace utopia {

    template<>
    class MultiReduce<TpetraVector, TRILINOS> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;
        static Scalar multi_min(const TpetraVector &t1, const TpetraVector &t2);
    };

}

#endif //UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP
