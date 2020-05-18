#include "utopia_kokkos_Eval_MultiReduce.hpp"

#include <utility>
#include "utopia_Tpetra_Vector.hpp"

#include "KokkosBlas1_dot.hpp"

namespace utopia {

    template <class Scalar, int N>
    class Tuple {
    public:
        Scalar data_[N];

        KOKKOS_INLINE_FUNCTION
        void init(const Scalar val) {
            for (int i = 0; i < N; i++) {
                data_[i] = val;
            }
        }

        template <class Op>
        KOKKOS_INLINE_FUNCTION void join(const Op &op, const Tuple &other) {
            for (int i = 0; i < N; ++i) {
                data_[i] = op.apply(data_[i], other.data_[i]);
            }
        }

        template <class Op>
        KOKKOS_INLINE_FUNCTION void join(const Op &op, const volatile Tuple &other) volatile {
            for (int i = 0; i < N; ++i) {
                auto v = data_[i];
                auto v2 = other.data_[i];
                data_[i] = op.apply(v, v2);
            }
        }

        KOKKOS_INLINE_FUNCTION
        Tuple operator+=(const Tuple &other) {
            for (int i = 0; i < N; ++i) {
                data_[i] = other.data_[i];
            }
        }

        KOKKOS_INLINE_FUNCTION
        volatile Tuple &operator+=(const volatile Tuple &other) volatile {
            for (int i = 0; i < N; ++i) {
                // auto v = data_[i];
                auto v2 = other.data_[i];
                data_[i] += v2;
            }

            return *this;
        }

        KOKKOS_INLINE_FUNCTION
        const Scalar &operator[](const int i) const { return data_[i]; }

        KOKKOS_INLINE_FUNCTION
        Scalar &operator[](const int i) { return data_[i]; }
    };

    template <class T, class Data, class Space, class Op, int N>
    struct MultiOpReducer {
    public:
        // Required
        using reducer = MultiOpReducer<T, Data, Space, Op, N>;
        typedef utopia::Tuple<T, N> TupleT;
        typedef Kokkos::View<TupleT *, Space, Kokkos::MemoryUnmanaged> result_view_type;

    private:
        Op op_;
        Tuple<Data, N> data_;
        T init_val_;

    public:
        KOKKOS_INLINE_FUNCTION
        MultiOpReducer(const Op &op, Tuple<Data, N> data, const T &init_val)
            : op_(op), data_(std::move(data)), init_val_(init_val) {}

        KOKKOS_INLINE_FUNCTION
        void join(volatile TupleT &val, const volatile TupleT &other) const {
            // Kokkos forces us to have the input values being declared volatile. Hence we need to make copies for the
            // reduction operations
            val.join(op_, other);
        }

        // Required
        KOKKOS_INLINE_FUNCTION
        void join(TupleT &dest, const TupleT &src) const { dest.join(op_, src); }

        KOKKOS_INLINE_FUNCTION void operator()(const int &i, TupleT &val) const {
            for (int valIdx = 0; valIdx < N; ++valIdx) {
                val[valIdx] = op_.apply(val[valIdx], data_[valIdx](i, 0));
            }
        }

        KOKKOS_INLINE_FUNCTION
        void init(TupleT &val) const { val.init(init_val_); }
    };

    template <class Vector, class Op>
    class KokkosEvalMultiReduce {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using ExecutionSpaceT = typename Vector::ExecutionSpace;

        inline static void eval(const Vector &t1,
                                const Vector &t2,
                                // const Op &,
                                const Scalar initial_value,
                                Tuple<Scalar, 2> &result) {
            using Data = decltype(t1.raw_type()->template getLocalView<ExecutionSpaceT>());

            Tuple<Data, 2> data;
            data[0] = t1.raw_type()->template getLocalView<ExecutionSpaceT>();
            data[1] = t2.raw_type()->template getLocalView<ExecutionSpaceT>();

            KokkosOp<Scalar, Op> kop;
            MultiOpReducer<Scalar, Data, ExecutionSpaceT, KokkosOp<Scalar, Op>, 2> reducer(kop, data, initial_value);

            Kokkos::parallel_reduce(data[0].extent(0), reducer, result);
        }
    };

    TpetraVector::Scalar MultiReduce<TpetraVector, TRILINOS>::multi_min(const TpetraVector &t1,
                                                                        const TpetraVector &t2) {
        using Scalar = Traits<TpetraVector>::Scalar;

        assert(t1.local_size() == t2.local_size());

        Scalar max_val = std::numeric_limits<Scalar>::max();

        Tuple<Scalar, 2> result{};

        KokkosEvalMultiReduce<TpetraVector, Min>::eval(t1, t2, max_val, result);

        const Scalar local_min = std::min(result[0], result[1]);
        return t1.comm().min(local_min);
    }

    //////////////////////////////////////////////////////////////////////////

    void EvalDots<TpetraVector, TRILINOS>::apply(const TpetraVector &v11,
                                                 const TpetraVector &v12,
                                                 Scalar &result1,
                                                 const TpetraVector &v21,
                                                 const TpetraVector &v22,
                                                 Scalar &result2) {
        using ExecutionSpaceT = TpetraVector::ExecutionSpace;
        using Data = decltype(v11.raw_type()->template getLocalView<ExecutionSpaceT>());

        const auto &comm = v11.comm();

        Data d11 = v11.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d12 = v12.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d21 = v21.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d22 = v22.raw_type()->template getLocalView<ExecutionSpaceT>();

        std::array<Scalar, 2> result{};
        if (v11.local_size() == v21.local_size()) {
            Tuple<Scalar, 2> tuple_result{};
            tuple_result.init(0.0);

            Kokkos::parallel_reduce(d11.extent(0),
                                    KOKKOS_LAMBDA(const int i, Tuple<Scalar, 2> &t) {
                                        t[0] += d11(i, 0) * d12(i, 0);
                                        t[1] += d21(i, 0) * d22(i, 0);
                                    },
                                    tuple_result);

            result[0] = tuple_result[0];
            result[1] = tuple_result[1];

        } else {
            result[0] = KokkosBlas::dot(Kokkos::subview(d11, Kokkos::ALL(), 0), Kokkos::subview(d12, Kokkos::ALL(), 0));

            result[1] = KokkosBlas::dot(Kokkos::subview(d21, Kokkos::ALL(), 0), Kokkos::subview(d22, Kokkos::ALL(), 0));
        }

        comm.sum(result);
        result1 = result[0];
        result2 = result[1];
    }

    void EvalDots<TpetraVector, TRILINOS>::apply(const TpetraVector &v11,
                                                 const TpetraVector &v12,
                                                 Scalar &result1,
                                                 const TpetraVector &v21,
                                                 const TpetraVector &v22,
                                                 Scalar &result2,
                                                 const TpetraVector &v31,
                                                 const TpetraVector &v32,
                                                 Scalar &result3) {
        using ExecutionSpaceT = TpetraVector::ExecutionSpace;
        using Data = decltype(v11.raw_type()->template getLocalView<ExecutionSpaceT>());

        const auto &comm = v11.comm();

        Data d11 = v11.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d12 = v12.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d21 = v21.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d22 = v22.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d31 = v31.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d32 = v32.raw_type()->template getLocalView<ExecutionSpaceT>();

        std::array<Scalar, 3> result{};

        if (v11.local_size() == v21.local_size() && v11.local_size() == v31.local_size()) {
            Tuple<Scalar, 3> tuple_result{};
            tuple_result.init(0.0);

            Kokkos::parallel_reduce(d11.extent(0),
                                    KOKKOS_LAMBDA(const int i, Tuple<Scalar, 3> &t) {
                                        t[0] += d11(i, 0) * d12(i, 0);
                                        t[1] += d21(i, 0) * d22(i, 0);
                                        t[2] += d31(i, 0) * d32(i, 0);
                                    },
                                    tuple_result);

            result[0] = tuple_result[0];
            result[1] = tuple_result[1];
            result[2] = tuple_result[2];

        } else {
            result[0] = KokkosBlas::dot(Kokkos::subview(d11, Kokkos::ALL(), 0), Kokkos::subview(d12, Kokkos::ALL(), 0));

            result[1] = KokkosBlas::dot(Kokkos::subview(d21, Kokkos::ALL(), 0), Kokkos::subview(d22, Kokkos::ALL(), 0));

            result[2] = KokkosBlas::dot(Kokkos::subview(d31, Kokkos::ALL(), 0), Kokkos::subview(d32, Kokkos::ALL(), 0));
        }

        comm.sum(result);
        result1 = result[0];
        result2 = result[1];
        result3 = result[2];
    }

}  // namespace utopia
