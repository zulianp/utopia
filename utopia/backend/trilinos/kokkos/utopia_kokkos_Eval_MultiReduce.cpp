#include "utopia_kokkos_Eval_MultiReduce.hpp"
#include "utopia_TpetraVector.hpp"

namespace utopia {

    template<typename Data, typename KokkosOp, typename Scalar, int N>
    struct MultiOpFunctor {
        KOKKOS_INLINE_FUNCTION void join(volatile Scalar &val[N], const volatile Scalar other[N]) const {
            // Kokkos forces us to have the input values being declared volatile. Hence we need to make copies for the reduction operations
            Scalar tmp1, tmp2;
            for(int d = 0; d < N; ++d) {
                tmp1 = val[d];
                tmp2 = other[d];
                val[d] = op_.apply(tmp1, tmp2);
            }
        }

        KOKKOS_INLINE_FUNCTION void operator()(const int& i, Scalar val[N]) const {
            for(int d = 0; d < N; ++d) {
                val[d] = op_.apply(val, data_[d](i, 0));
            }
        }

        KOKKOS_INLINE_FUNCTION void init(Scalar val[N]) const
        {
            for(int i = 0; i < N; ++i) {
                val[i] = initial_value_[i];
            }
        }

        MultiOpFunctor(
            const KokkosOp &op,
            const Data &data[N],
            const Scalar &initial_value[N])
        : op_(op)
        {
            for(int i = 0; i < N; ++i) {
                data_[i]          = data[i];
                initial_value_[i] = initial_value[i];
            }
        }

    private:
        const KokkosOp op_;
        const Data data_[N];
        const Scalar initial_value_[N];
    };

    template<class Vector, class Op>
    class KokkosEvalMultiReduce {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        template<int N>
        inline static void eval(
            const Vector &vec,
            const Op op,
            const Scalar &initial_value[N],
            Scalar &result[N],
            )
        {
            using ExecutionSpaceT = typename Vector::vector_type::execution_space;
            using Scalar = typename Vector::Scalar;
            using Data = decltype(vec.raw_type()->template getLocalView<ExecutionSpaceT>());

            assert(!vec.empty());
            auto data = vec.raw_type()->template getLocalView<ExecutionSpaceT>();


            KokkosOp<Scalar, Op> kop;
            MultiOpFunctor<Data, KokkosOp<Scalar, Op>, Scalar, N> functor{ kop, data, initial_value };
            Kokkos::parallel_reduce(data.extent(0), functor, result);

            return ret;
        }
    };


    TpetraVector::Scalar MultiReduce<TpetraVector, TRILINOS>::multi_min(
        const TpetraVector &t1,
        const TpetraVector &t2
    )
    {

    }

}
