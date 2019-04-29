#ifndef UTOPIA_KOKKOS_EVAL_DISTANCE_HPP
#define UTOPIA_KOKKOS_EVAL_DISTANCE_HPP

#include "utopia_kokkos_Operations.hpp"
#include "utopia_Norm.hpp"

#include <Kokkos_Core.hpp>
#include <iostream>
#include <cassert>

namespace utopia {

    template<class Vector, int NormType>
    class KokkosEvalDistance {};

    template<class Vector>
    class KokkosEvalDistance<Vector, 2> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        inline static Scalar finalize(const Scalar &value)
        {
            return std::sqrt(value);
        }

        inline static Scalar apply(const Vector &left, const Vector &right, const bool finalize_reduction = true)
        {
            using ExecutionSpaceT = typename Vector::vector_type::execution_space;
            using Scalar = typename Vector::Scalar;

            assert(!left.empty());
            assert(left.size() == right.size());
            assert(left.local_size() == right.local_size());


            auto l_data = left.implementation().template getLocalView<ExecutionSpaceT>();
            auto r_data = right.implementation().template getLocalView<ExecutionSpaceT>();

            Scalar ret = 0.;

            Kokkos::parallel_reduce(l_data.extent(0), KOKKOS_LAMBDA(const int i, Scalar &val) {
                const Scalar x = l_data(i, 0) - r_data(i, 0);
                val += x * x;
            }, ret);

            if(finalize_reduction) {
                return finalize(ret);
            } else {
                return ret;
            }
        }
    };

    template<class Vector>
    class KokkosEvalDistance<Vector, 1> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        inline static Scalar finalize(const Scalar &value)
        {
            return value;
        }

        inline static Scalar apply(const Vector &left, const Vector &right, const bool)
        {
            using ExecutionSpaceT = typename Vector::vector_type::execution_space;
            using Scalar = typename Vector::Scalar;

            assert(!left.empty());
            assert(left.size() == right.size());
            assert(left.local_size() == right.local_size());


            auto l_data = left.implementation().template getLocalView<ExecutionSpaceT>();
            auto r_data = right.implementation().template getLocalView<ExecutionSpaceT>();

            Scalar ret = 0.;

            Kokkos::parallel_reduce(l_data.extent(0), KOKKOS_LAMBDA(const int i, Scalar &val) {
                const Scalar x = l_data(i, 0) - r_data(i, 0);
                val += Kokkos::Details::ArithTraits<Scalar>::abs(x);
            }, ret);

            return ret;
        }
    };

    template<class Vector>
    class KokkosEvalDistance<Vector, INFINITY_NORM_TAG> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        inline static Scalar finalize(const Scalar &value)
        {
            return value;
        }

        inline static Scalar apply(const Vector &left, const Vector &right, const bool)
        {
            using ExecutionSpaceT = typename Vector::vector_type::execution_space;
            using Scalar = typename Vector::Scalar;

            assert(!left.empty());
            assert(left.size() == right.size());
            assert(left.local_size() == right.local_size());


            auto l_data = left.implementation().template getLocalView<ExecutionSpaceT>();
            auto r_data = right.implementation().template getLocalView<ExecutionSpaceT>();

            Scalar ret = 0.;

            Kokkos::parallel_reduce(l_data.extent(0), KOKKOS_LAMBDA(const int i, Scalar &val) {
                const Scalar x = l_data(i, 0) - r_data(i, 0);
                val = KokkosOp<Scalar, Max>::apply(val, Kokkos::Details::ArithTraits<Scalar>::abs(x));
            }, ret);

            return ret;
        }
    };

}

#endif //UTOPIA_KOKKOS_EVAL_DISTANCE_HPP
