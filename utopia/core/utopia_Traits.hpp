//
// Created by Patrick Zulian on 19/05/15.
//

#ifndef utopia_utopia_TRAITS_HPP
#define utopia_utopia_TRAITS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include <vector>



namespace utopia {
    static const int HOMEMADE = -1;
    static const int BLAS = 1;
    static const int PETSC = 10;
    static const int CUDA = 100;
    static const int OPENCL_TAG = 1000;

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


    template<class Traits, int Order, int Sparsity = FillType::DENSE>
    class TensorQuery {};

    template<class Traits, int Sparsity>
    class TensorQuery<Traits, 0, Sparsity> {
    public:
        typedef typename Traits::Scalar Type;
    };

    template<class Traits, int Sparsity>
    class TensorQuery<Traits, 1, Sparsity> {
    public:
        typedef typename Traits::Vector Type;
    };

    template<class Traits>
    class TensorQuery<Traits, 2, FillType::DENSE> {
    public:
        typedef typename Traits::Matrix Type;
    };

    template<class Traits>
    class TensorQuery<Traits, 2, FillType::SPARSE> {
    public:
        typedef typename Traits::SparseMatrix Type;
    };

    template<class Left, class Right>
    class MostDescriptive {
        public:
        typedef Left Type;
    };


    template<class Left, class Right, class Default,
             int SparsityLeft  = utopia::Traits<Left>::FILL_TYPE, 
             int SparsityRight = utopia::Traits<Right>::FILL_TYPE>
    class ChooseType {
    public:
        typedef Default Type;
    };

    template<class Left, class Right, class Default, int SparsityRight>
    class ChooseType<Left, Right, Default, FillType::DELEGATE, SparsityRight> {
    public:
        typedef Right Type;
    };

    template<class Left, class Right, class Default, int SparsityLeft>
    class ChooseType<Left, Right, Default, SparsityLeft,  FillType::DELEGATE> {
    public:
        typedef Left Type;
    };

    template<class Left, class Right, class Default>
    class ChooseType<Left, Right, Default, FillType::DELEGATE,  FillType::DELEGATE> {
    public:
        typedef Default Type;
    };


    template<class Left, class Right>
    class MostDescriptive<Left, Number<Right> > {
    public:
        typedef Left Type;
    };

    template<class Left, class Right>
    class MostDescriptive<Number<Left>, Right > {
    public:
        typedef Right Type;
    };
    
    

    // #define WRAPPER_TYPE(Traits, Expr) utopia::Wrapper<typename utopia::TensorQuery<Traits, Expr::Order, FillTypeQuery<Expr>::value>::Type, Expr::Order>

    /*!
     * @class Traits
     */
    template<class T>
    class Traits : public Traits<typename T::Implementation> {
    public:
        //You need to implement this for your type. See BLASTraits as an example
        //typedef double Scalar;
        // enum {
        //     FILL_TYPE = FillType::DENSE
        // };

        enum {
            Order = T::Order
        };
    };

    template<>
    class Traits<double>  {
    public:
        typedef double Scalar;

        enum {
            Order = 0
        };

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    template<>
    class Traits<float>  {
    public:
        typedef float Scalar;

        enum {
            Order = 0
        };

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    template<>
    class Traits<int>  {
    public:
        typedef int Scalar;

        enum {
            Order = 0
        };

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };


    template<class Tensor, class TraitsT = Traits<Tensor> >
    struct is_sparse {
        enum {
            value = TraitsT::FILL_TYPE == FillType::SPARSE
        };
    };


    template<class Traits, class Expr, int Sparsity = utopia::Traits<Expr>::FILL_TYPE>
    class TypeAndFill {
    public:
        typedef typename TensorQuery<Traits, Expr::Order, Sparsity>::Type Type;
        // typedef typename TensorQuery<Traits, Expr::Order>::Type Type;
    };

    #define EXPR_TYPE(ImplementationTraits, Expr) typename utopia::TypeAndFill<ImplementationTraits, Expr>::Type

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

#define UTOPIA_MAKE_TRAITS(TensorType, TraitsType)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultDelegateTraits {}; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultDelegateTraits {}; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultDelegateTraits {}

#define UTOPIA_MAKE_TRAITS_SPARSE(TensorType, TraitsType)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultSparseTraits { }; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultSparseTraits { }; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultSparseTraits { }

#define UTOPIA_MAKE_TRAITS_DENSE(TensorType, TraitsType)  \
    template<> class Traits<TensorType> : public TraitsType, public DefaultDenseTraits { }; \
    template<> class Traits<const TensorType &> : public TraitsType, public DefaultDenseTraits { }; \
    template<> class Traits<TensorType &> : public TraitsType, public DefaultDenseTraits { } 

#define UTOPIA_MAKE_TRAITS_TPL_1(TensorType, TraitsType)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultDelegateTraits {}; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultDelegateTraits {}; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultDelegateTraits {}

#define UTOPIA_MAKE_TRAITS_SPARSE_TPL_1(TensorType, TraitsType)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultSparseTraits { }; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultSparseTraits { }; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultSparseTraits { }   

#define UTOPIA_MAKE_TRAITS_DENSE_TPL_1(TensorType, TraitsType)  \
    template<typename T> class Traits< TensorType<T> > : public TraitsType<T>, public DefaultDenseTraits { }; \
    template<typename T> class Traits<const TensorType<T> &> : public TraitsType<T>, public DefaultDenseTraits { }; \
    template<typename T> class Traits< TensorType<T> &> : public TraitsType<T>, public DefaultDenseTraits {  }     

    //Always to be called within the utopia namespace
#define UTOPIA_MAKE_PARALLEL_TRAITS(TensorType)  template<> class is_parallel<TensorType> { public: enum { value = true}; }
}

#endif //utopia_utopia_TRAITS_HPP
