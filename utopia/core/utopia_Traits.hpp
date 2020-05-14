//
// Created by Patrick Zulian on 19/05/15.
//

#ifndef utopia_utopia_TRAITS_HPP
#define utopia_utopia_TRAITS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include <vector>



namespace utopia {
    static const int INVALID_BACKEND = -100;
    static const int HOMEMADE = -1;
    static const int BLAS = 1;
    static const int PETSC = 10;
    static const int CUDA = 100;
    static const int OPENCL_TAG = 1000;
    static const int TRILINOS = 10000;
    static const int KOKKOS   = 100000;
    static const int SERIAL_HOMEMADE = 200000;
    static const int PETSC_EXPERIMENTAL = -1000;

    class FillType {
    public:
        /////////////////////////////////////////////////
        // static const int ZERO               = 0;

        /////////////////////////////////////////////////
        static const int SPARSE             = 3;
        //for the moment the next ones are subclasses of SPARSE
        // static const int DIAGONAL           = 1;
        // static const int TRI_DIAGONAL       = 2;
        /////////////////////////////////////////////////

        /////////////////////////////////////////////////
        static const int DENSE              = 12;

        // For the moment the next ones are subclasses of DENSE
        // static const int LOWER_TRIANGULAR   = 4;
        // static const int UPPER_TRIANGULAR   = 8;
        static const int SCALAR             = 16;
        static const int DELEGATE           = 32;
        static const int POLYMORPHIC        = 64;

        /////////////////////////////////////////////////
    };

    template<int Type>
    class Fill2String {
    public:
        constexpr static const char * Value()
        {
            return "Undefined";
        }
    };

    template<>
    class Fill2String<FillType::SPARSE> {
    public:
        constexpr static const char * Value()
        {
            return "Sparse";
        }
    };

    template<>
    class Fill2String<FillType::DENSE> {
    public:
        constexpr static const char * Value()
        {
            return "Dense";
        }
    };

    template<>
    class Fill2String<FillType::DELEGATE> {
    public:
        constexpr static const char * Value()
        {
            return "";
        }
    };

    template<>
    class Fill2String<FillType::SCALAR> {
    public:
        constexpr static const char * Value()
        {
            return "Scalar";
        }
    };


    template<>
    class Fill2String<FillType::POLYMORPHIC> {
    public:
        constexpr static const char * Value()
        {
            return "Polymorphic";
        }
    };

    template<typename T>
    class Unwrap {
    public:
        using Type = T;
    };

    template<class Derived, int Order>
    class Unwrap<Tensor<Derived, Order>> {
    public:
        using Type = Derived;
    };

    template<class Traits, int Order, int Sparsity = FillType::DENSE>
    class TensorQuery {};

    template<class Traits, int Sparsity>
    class TensorQuery<Traits, 0, Sparsity> {
    public:
        using Type = utopia::Number<typename Traits::Scalar>;
    };

    template<class Traits, int Sparsity>
    class TensorQuery<Traits, 1, Sparsity> {
    public:
        using Type = typename Traits::Vector;
    };

    template<class Traits>
    class TensorQuery<Traits, 2, FillType::DENSE> {
    public:
        using Type = typename Traits::Matrix;
    };

    template<class Traits>
    class TensorQuery<Traits, 2, FillType::SPARSE> {
    public:
        using Type = typename Traits::SparseMatrix;
    };

    template<class Traits>
    class TensorQuery<Traits, 2, FillType::POLYMORPHIC> {
    public:
        using Type = typename Traits::PolymorphicMatrix;
    };


    template<class Traits, int Sparsity>
    class TensorQuery<Traits, 4, Sparsity> {
    public:
        //FIXME dummy for 4th order tensor
        typedef typename TensorQuery<Traits, 2, Sparsity>::Type Type;
    };


    template<class Left, class Right>
    class MostDescriptive {
    public:
        using Type = typename Unwrap<Left>::Type;
    };


    template<class Left, class Right, class Default,
             int SparsityLeft  = utopia::Traits<Left>::FILL_TYPE,
             int SparsityRight = utopia::Traits<Right>::FILL_TYPE>
    class ChooseType {
    public:
        using Type = typename Unwrap<Default>::Type;
    };

    template<class Left, class Right, class Default, int SparsityRight>
    class ChooseType<Left, Right, Default, FillType::DELEGATE, SparsityRight> {
    public:
        using Type = typename Unwrap<Right>::Type;
    };

    template<class Left, class Right, class Default, int SparsityLeft>
    class ChooseType<Left, Right, Default, SparsityLeft,  FillType::DELEGATE> {
    public:
        using Type = typename Unwrap<Left>::Type;
    };

    template<class Left, class Right, class Default>
    class ChooseType<Left, Right, Default, FillType::DELEGATE,  FillType::DELEGATE> {
    public:
        using Type = typename Unwrap<Default>::Type;
    };


    template<class Left, class Right>
    class MostDescriptive<Left, Number<Right> > {
    public:
        using Type = typename Unwrap<Left>::Type;
    };

    template<class Left, class Right>
    class MostDescriptive<Number<Left>, Right > {
    public:
        using Type = typename Unwrap<Right>::Type;
    };

    template<typename Left, typename Right>
    class MostDescriptive<Number<Left>, Number<Right> > {
    public:
        // typedef decltype(Left(0) + Right(0)) Type;
        using Type = utopia::Number<decltype(Left(0) + Right(0))>;
    };


    // #define WRAPPER_TYPE(Traits, Expr) utopia::Wrapper<typename utopia::TensorQuery<Traits, Expr::Order, FillTypeQuery<Expr>::value>::Type, Expr::Order>

    /*!
     * @class Traits
     */
    // template<class T>
    // class Traits : public Traits<typename T::Implementation> {
    // public:
    //     //You need to implement this for your type. See BLASTraits as an example
    //     //typedef double Scalar;
    //     // enum {
    //     //     FILL_TYPE = FillType::DENSE
    //     // };

    //     static const int Order = T::Order;
    // };

    /*!
     * @class Traits
     */
    template<class T>
    class Traits<const T> : public Traits<T> {
    public:
    };

    // template<typename T>
    // class Traits< Number<T> > : public Traits<T> {};


    template<>
    class Traits<double>  {
    public:
        using Scalar = double;

        static const int Backend = INVALID_BACKEND;
        static const int Order = 0;

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    template<>
    class Traits<float>  {
    public:
        using Scalar = float;

        static const int Backend = INVALID_BACKEND;
        static const int Order = 0;

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    template<>
    class Traits<int>  {
    public:
        using Scalar = int;

        static const int Backend = INVALID_BACKEND;
        static const int Order = 0;

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };


    template<class Tensor, class TraitsT = Traits<Tensor> >
    struct is_sparse {
        enum {
            value = TraitsT::FILL_TYPE == FillType::SPARSE || TraitsT::FILL_TYPE == FillType::POLYMORPHIC
        };
    };

    template<class Tensor, class TraitsT = Traits<Tensor> >
    struct is_polymorhic {
        enum {
            value = TraitsT::FILL_TYPE == FillType::POLYMORPHIC
        };
    };

    template<class Tensor, class TraitsT = Traits<Tensor> >
    struct is_dense {
        enum {
            value = TraitsT::FILL_TYPE == FillType::DENSE
        };
    };

    template<class Tensor, class TraitsT = Traits<Tensor> >
    struct is_dense_or_polymorphic {
        enum {
            value = is_dense<Tensor>::value || is_polymorhic<Tensor>::value
        };
    };


    template<class Traits, class Expr, int Sparsity = utopia::Traits<Expr>::FILL_TYPE>
    class TypeAndFill {
    public:
        typedef typename TensorQuery<Traits, Expr::Order, Sparsity>::Type Type;
        // typedef typename TensorQuery<Traits, Expr::Order>::Type Type;
    };

    #define EXPR_TYPE(ImplementationTraits, Expr) typename utopia::TypeAndFill<ImplementationTraits, Expr>::Type

    template<int Order_>
    class DefaultTensorTraits {
    public:
        static const int Order = Order_;
    };

    class DefaultDenseTraits {
    public:
        enum {
            FILL_TYPE = FillType::DENSE
        };
    };

    class DefaultDelegateTraits {
    public:
        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    class DefaultSparseTraits {
    public:
        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    class DefaultPolymorphicTraits {
    public:
        enum {
            FILL_TYPE = FillType::POLYMORPHIC
        };
    };

#define UTOPIA_MAKE_TRAITS(TensorType, TraitsType, Order)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}

#define UTOPIA_MAKE_TRAITS_SPARSE(TensorType, TraitsType, Order)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }

#define UTOPIA_MAKE_TRAITS_DENSE(TensorType, TraitsType, Order)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultDenseTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultDenseTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultDenseTraits, public DefaultTensorTraits<Order> { }

#define UTOPIA_MAKE_TRAITS_POLYMORPHIC(TensorType, TraitsType, Order)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultPolymorphicTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultPolymorphicTraits, public DefaultTensorTraits<Order> { }; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultPolymorphicTraits, public DefaultTensorTraits<Order> { }

#define UTOPIA_MAKE_TRAITS_TPL_1(TensorType, TraitsType, Order)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultDelegateTraits, public DefaultTensorTraits<Order> {}

#define UTOPIA_MAKE_TRAITS_SPARSE_TPL_1(TensorType, TraitsType, Order)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultSparseTraits, public DefaultTensorTraits<Order> { }

#define UTOPIA_MAKE_TRAITS_DENSE_TPL_1(TensorType, TraitsType, Order)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultDenseTraits, public DefaultTensorTraits<Order> { }; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultDenseTraits, public DefaultTensorTraits<Order> { }; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultDenseTraits, public DefaultTensorTraits<Order> {  }

    //Always to be called within the utopia namespace
#define UTOPIA_MAKE_PARALLEL_TRAITS(TensorType)  template<> class is_parallel<TensorType> { public: enum { value = true}; }
}

#endif //utopia_utopia_TRAITS_HPP
