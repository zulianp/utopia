#ifndef UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP
#define UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP

#include "utopia_trilinos_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_MultiReduce.hpp"
#include "utopia_DotVecVecs.hpp"

namespace utopia {

    template<>
    class MultiReduce<TpetraVector, TRILINOS> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;

        static Scalar multi_min(const TpetraVector &t1, const TpetraVector &t2);
    };


    template<>
    class EvalDots<TpetraVector, TRILINOS> : public EvalDots<TpetraVector, HOMEMADE> {
    public:
        using Scalar = typename utopia::Traits<TpetraVector>::Scalar;

        // static void apply(
        //     const TpetraVector &v1,
        //     const std::vector<std::shared_ptr<TpetraVector> > &vectors,
        //     std::vector<Scalar> & results
        // );

        // static void apply(
        //     const TpetraVector &v1,
        //     const std::vector<TpetraVector> &vectors,
        //     std::vector<Scalar> & results
        // );

        static void apply(
            const TpetraVector &v11,
            const TpetraVector &v12,
            Scalar & result1,
            const TpetraVector &v21,
            const TpetraVector &v22,
            Scalar & result2
        );

        static void apply(
            const TpetraVector &v11,
            const TpetraVector &v12,
            Scalar & result1,
            const TpetraVector &v21,
            const TpetraVector &v22,
            Scalar & result2,
            const TpetraVector &v31,
            const TpetraVector &v32,
            Scalar & result3
        );
    };

}

#endif //UTOPIA_KOKKOS_EVAL_MULTI_REDUCE_HPP
