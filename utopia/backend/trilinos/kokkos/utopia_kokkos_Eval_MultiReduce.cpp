#include "utopia_kokkos_Eval_MultiReduce.hpp"
#include "utopia_Tpetra_Vector.hpp"

#include "KokkosBlas1_dot.hpp"

namespace utopia {

    template<class Scalar, int N>
    class Tuple {
    public:
      Scalar data_[N];

      KOKKOS_INLINE_FUNCTION
      Tuple() {
         // init();
      }

      KOKKOS_INLINE_FUNCTION
      Tuple(const Tuple & rhs) {
         for (int i = 0; i < N; i++ ){
            data_[i] = rhs.data_[i];
         }
      }

      KOKKOS_INLINE_FUNCTION  // initialize data_ to 0
      void init(const Scalar val) {
          for (int i = 0; i < N; i++ ) { data_[i] = val; }
       }

   template<class Op>
    void join(const Op &op, const Tuple &other)
    {
        for(int i = 0; i < N; ++i) {
            data_[i] = op.apply(data_[i], other.data_[i]);
        }
    }

    template<class Op>
     void join(const Op &op, const volatile Tuple &other) volatile
     {
         for(int i = 0; i < N; ++i) {
             data_[i] = op.apply(data_[i], other.data_[i]);
         }
     }

    inline const Scalar &operator[](const int i) const
    {
        // assert(i < N);
        return data_[i];
    }

    inline Scalar &operator[](const int i)
    {
        // assert(i < N);
        return data_[i];
    }

      // KOKKOS_INLINE_FUNCTION
      // Tuple& operator += (const Tuple& src) {
      //   for ( int i = 0; i < N; i++ ) {
      //      data_[i]+=src.data_[i];
      //   }
      //   return *this;
      // }

      // KOKKOS_INLINE_FUNCTION
      // void operator += (const volatile Tuple& src) volatile {
      //   for ( int i = 0; i < N; i++ ) {
      //     data_[i]+=src.data_[i];
      //   }
      // }

    };

    template<class T, class Space, class Op, int N>
    struct MultiOpReducer {
    public:
      //Required
      typedef MultiOpReducer reducer;
      typedef utopia::Tuple<T,N> value_type;
      typedef Kokkos::View<value_type*, Space, Kokkos::MemoryUnmanaged> result_view_type;

    private:
      Op op_;
      T init_val_;
      value_type & value;
      
    public:

      KOKKOS_INLINE_FUNCTION
      MultiOpReducer(
        const Op &op,
        const T &init_val,
        value_type& value_): op_(op), init_val_(init_val), value(value_) {}

      //Required
      KOKKOS_INLINE_FUNCTION
      void join(value_type& dest, const value_type& src)  const {
        dest.join(op_, src);
      }

      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type& dest, const volatile value_type& src) const {
        dest += src;
      }

      KOKKOS_INLINE_FUNCTION
      void init( value_type& val)  const {
        val.init(init_val_);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& reference() const {
        return value;
      }

      KOKKOS_INLINE_FUNCTION
      result_view_type view() const {
        return result_view_type(&value);
      }

      KOKKOS_INLINE_FUNCTION
      bool references_scalar() const {
        return true;
      }
    };


    template<class Vector, class Op>
    class KokkosEvalMultiReduce {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using ExecutionSpaceT = typename Vector::vector_type::execution_space;

        inline static void eval(
            const Vector &t1,
            const Vector &t2,
            // const Op &,
            const Scalar initial_value,
            Tuple<Scalar, 2> &result
            )
        {
            using Data = decltype(t1.raw_type()->template getLocalView<ExecutionSpaceT>());

            Data d1 = t1.raw_type()->template getLocalView<ExecutionSpaceT>();
            Data d2 = t2.raw_type()->template getLocalView<ExecutionSpaceT>();

            KokkosOp<Scalar, Op> kop;
            MultiOpReducer<Scalar, ExecutionSpaceT, KokkosOp<Scalar, Op>, 2> reducer(kop, initial_value, result);
            
            Kokkos::parallel_reduce(d1.extent(0),
                KOKKOS_LAMBDA(const int i, Tuple<Scalar, 2> &tuple) {
                    tuple[0] = d1(i, 0);
                    tuple[1] = d2(i, 0);
                }, reducer);
        }
    };

    TpetraVector::Scalar MultiReduce<TpetraVector, TRILINOS>::multi_min(
        const TpetraVector &t1,
        const TpetraVector &t2
    )
    {
        using Scalar = Traits<TpetraVector>::Scalar;

        assert(t1.local_size() == t2.local_size());

        Scalar max_val = std::numeric_limits<Scalar>::max();

        Tuple<Scalar, 2> result;
       
        KokkosEvalMultiReduce<TpetraVector, Min>::eval(
            t1,
            t2,
            max_val,
            result
        );

        const Scalar local_min = std::min(result[0], result[1]);
        return t1.comm().min(local_min);
    }

    //////////////////////////////////////////////////////////////////////////

    // void EvalDots<TpetraVector, TRILINOS>::apply(
    //     const TpetraVector &v1,
    //     const std::vector<std::shared_ptr<TpetraVector> > &vectors,
    //     std::vector<Scalar> & results)
    // {
    //     typename utopia::Traits<Vector>::SizeType n =  vectors.size();

    //     if(n!=static_cast<SizeType>(results.size()))
    //         results.resize(n);

    //     std::vector<Vec> vecs(n);

    //     for(auto i=0; i < static_cast<SizeType>(vectors.size()); i++)
    //         vecs[i]=(raw_type(*vectors[i]));

    //      VecMDot(raw_type(v1), n, vecs.data(), results.data());
    // }

    // void EvalDots<TpetraVector, TRILINOS>::apply(
    //     const TpetraVector &v1,
    //     const std::vector<TpetraVector> &vectors,
    //     std::vector<Scalar> & results)
    // {
    //     std::vector<Vec> vecs;
    //     for(auto i=0; i < static_cast<SizeType>(vectors.size()); i++)
    //     {
    //         if(!empty(vectors[i]))
    //             vecs.push_back(raw_type(vectors[i]));
    //     }

    //     typename utopia::Traits<Vector>::SizeType n =  vecs.size();

    //     if(n != vectors.size() || n != results.size())
    //     {
    //         std::vector<Scalar> result_new(n);
    //         VecMDot(raw_type(v1), n, vecs.data(), result_new.data());

    //         for(auto i=0; i < n; i++)
    //             results[i] = result_new[i];
    //     }
    //     else
    //     {
    //         if(n!=results.size())
    //             results.resize(n);

    //         VecMDot(raw_type(v1), n, vecs.data(), results.data());
    //     }
    // }

    void EvalDots<TpetraVector, TRILINOS>::apply(
        const TpetraVector &v11,
        const TpetraVector &v12,
        Scalar &result1,
        const TpetraVector &v21,
        const TpetraVector &v22,
        Scalar &result2)
    {

        using ExecutionSpaceT = typename TpetraVector::vector_type::execution_space;
        using Data = decltype(v11.raw_type()->template getLocalView<ExecutionSpaceT>());

        const auto &comm = v11.comm();

        Data d11 = v11.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d12 = v12.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d21 = v21.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d22 = v22.raw_type()->template getLocalView<ExecutionSpaceT>();

        std::array<Scalar, 2> result;
        result[0] = KokkosBlas::dot(
                Kokkos::subview(d11, Kokkos::ALL(), 0),
                Kokkos::subview(d12, Kokkos::ALL(), 0)
        );

        result[1] = KokkosBlas::dot(
                Kokkos::subview(d21, Kokkos::ALL(), 0),
                Kokkos::subview(d22, Kokkos::ALL(), 0)
        );

        comm.sum(result);
        result1 = result[0];
        result2 = result[1];
    }

    void EvalDots<TpetraVector, TRILINOS>::apply(
        const TpetraVector &v11,
        const TpetraVector &v12,
        Scalar &result1,
        const TpetraVector &v21,
        const TpetraVector &v22,
        Scalar &result2,
        const TpetraVector &v31,
        const TpetraVector &v32,
        Scalar &result3)
    {
        using ExecutionSpaceT = typename TpetraVector::vector_type::execution_space;
        using Data = decltype(v11.raw_type()->template getLocalView<ExecutionSpaceT>());

        const auto &comm = v11.comm();

        Data d11 = v11.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d12 = v12.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d21 = v21.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d22 = v22.raw_type()->template getLocalView<ExecutionSpaceT>();

        Data d31 = v31.raw_type()->template getLocalView<ExecutionSpaceT>();
        Data d32 = v32.raw_type()->template getLocalView<ExecutionSpaceT>();

        std::array<Scalar, 3> result;
        result[0] = KokkosBlas::dot(
                Kokkos::subview(d11, Kokkos::ALL(), 0),
                Kokkos::subview(d12, Kokkos::ALL(), 0)
        );

        result[1] = KokkosBlas::dot(
                Kokkos::subview(d21, Kokkos::ALL(), 0),
                Kokkos::subview(d22, Kokkos::ALL(), 0)
        );

        result[2] = KokkosBlas::dot(
                Kokkos::subview(d31, Kokkos::ALL(), 0),
                Kokkos::subview(d32, Kokkos::ALL(), 0)
        );

        comm.sum(result);
        result1 = result[0];
        result2 = result[1];
        result3 = result[1];
    }

}
