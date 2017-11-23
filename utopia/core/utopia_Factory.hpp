//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_FACTORY_HPP
#define UTOPIA_UTOPIA_FACTORY_HPP

#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"

namespace utopia {

    template<class Type>
    class FactoryTraits {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "TODO";
        }

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    class Identity {};

    template<>
    class FactoryTraits<Identity> {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "Identity";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    class LocalIdentity {};
    class Zeros {};
    class LocalZeros {};



    template<>
    class FactoryTraits<LocalIdentity> {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "LocalIdentity";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };


    template<>
    class FactoryTraits<LocalZeros> {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "LocalZeros";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };


    template<>
    class FactoryTraits<Zeros> {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "Zeros";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    template<typename T>
    class Values {
    public:

        Values()  {};
        Values(T value) : _value(value) {};
        template<typename OtherT>
        Values(const Values<OtherT> & other) {_value = other.value();}

        inline const T &value() const
        {
            return _value;
        }
        inline void setValue(T value)
        {
            _value = value;
        }

    private:
        T _value;
    };

    template<typename _Scalar>
    class FactoryTraits< Values<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static constexpr const char * getClass()
        {
            return "Values";
        }

        enum {
            FILL_TYPE = FillType::DENSE
        };
    };


    template<typename T>
    class LocalValues {
    public:

        LocalValues()  {};
        LocalValues(T value) : _value(value) {};
        template<typename OtherT>
        LocalValues(const LocalValues<OtherT> & other) {_value = other.value();}

        inline const T &value() const
        {
            return _value;
        }
        inline void setValue(T value)
        {
            _value = value;
        }

    private:
        T _value;
    };

    template<typename _Scalar>
    class FactoryTraits< LocalValues<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static constexpr const char * getClass()
        {
            return "LocalValues";
        }

        enum {
            FILL_TYPE = FillType::DENSE
        };
    };


    template<typename T>
    class NNZ {
    public:
    
        NNZ()  {};
        NNZ(T nnz) : _nnz(nnz) {};
        template<typename OtherT>
        NNZ(const NNZ<OtherT>& other) {_nnz = other.nnz();}

        inline const T &nnz() const
        {
            return _nnz;
        }
        inline void setSparse(T nnz)
        {
            _nnz = nnz;
        }

    private:
        T _nnz;
    };

    template<typename _Scalar>
    class FactoryTraits< NNZ<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static constexpr const char * getClass()
        {
            return "NNZ";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };


    template<typename T>
    class LocalNNZ {
    public:
  

        LocalNNZ()  {};
        LocalNNZ(T nnz) : _nnz(nnz) {};
        template<typename OtherT>
        LocalNNZ(const LocalNNZ<OtherT>& other) {_nnz = other.nnz();}

        inline const T &nnz() const
        {
            return _nnz;
        }
        inline void setSparse(T nnz)
        {
            _nnz = nnz;
        }

    private:
        T _nnz;
    };


    template<typename _Scalar>
    class FactoryTraits< LocalNNZ<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static constexpr const char * getClass()
        {
            return "LocalNNZ";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    template<typename T>
    class LocalRowNNZ {
    public:
     
        LocalRowNNZ()  {};
        LocalRowNNZ(T nnz) : _nnz(nnz) {};
        template<typename OtherT>
        LocalRowNNZ(const LocalRowNNZ<OtherT>& other) {_nnz = other.nnz();}

        inline const T &nnz() const
        {
            return _nnz;
        }
        inline void setSparse(T nnz)
        {
            _nnz = nnz;
        }

    private:
        T _nnz;
    };

    
    template<typename _Scalar>
    class FactoryTraits< LocalRowNNZ<_Scalar> > {
    public:
        typedef _Scalar Scalar;

        static constexpr const char * getClass()
        {
            return "LocalRowNNZ";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };


    class Resize { };


    template<>
    class FactoryTraits<Resize> {
    public:
        typedef double Scalar;

        static constexpr const char * getClass()
        {
            return "Resize";
        }

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };


    template<class Type, int _Order>
    class Factory : public Expression< Factory<Type, _Order> > {
    public:
        static const int Order = _Order;

        enum {
            StoreAs = UTOPIA_BY_VALUE
        };

        typedef typename FactoryTraits<Type>::Scalar Scalar;

        inline const Size &size() const
        {
            return _size;
        }

        inline const Type &type() const
        {
            return _type;
        }
        Factory(const Size &size, const Type type = Type())
                : _size(size), _type(type)
        {}

        inline std::string getClass() const
        {
            return "Factory(" + std::string(FactoryTraits<Type>::getClass()) + ")";
        }

        virtual ~Factory() {}


    private:
        Size _size;
        Type _type;
    };


    template<class Type, int Order>
    class Traits< Factory<Type, Order> > {
    public:
        typedef typename utopia::FactoryTraits<Type>::Scalar Scalar;

        enum {
            FILL_TYPE = FactoryTraits<Type>::FILL_TYPE
        };
    };

    template<class Type, int Order>
    inline const Size &size(const Factory<Type, Order> &expr)
    {
        return expr.size();
    }


    /**    @defgroup factory Factory
     *      @brief  Factory methods allow for creating basic tensor in an easy way
     */


    /** \addtogroup Global
     * @brief Manipulation with objects on global adress space.
     * @ingroup factory
     *  @{
     */

    /// Returns identity matrix  \f$ I^{row \times cols}  \f$. 
    inline Factory<Identity, 2> identity(const Size::SizeType rows, const Size::SizeType cols)
    {
        return Factory<Identity, 2>(Size({rows, cols}));
    }
    /// Returns identity matrix  \f$ I^{size_0 \times size_1}  \f$. 
    inline Factory<Identity, 2> identity(const Size &size)
    {
        return Factory<Identity, 2>(size);
    }

    /// Returns global \f$ 0^{rows \times rows}  \f$. 
    inline Factory<Zeros, 2> zeros(const Size::SizeType rows, const Size::SizeType cols)
    {
        return Factory<Zeros, 2>(Size({rows, cols}));
    }

    /// Returns global \f$ 0^{n \times 1}  \f$. 
    inline Factory<Zeros, 1> zeros(const Size::SizeType n)
    {
        return Factory<Zeros, 1>(Size({n}));
    }

    /// Returns global \f$ 0^{size \times size}  \f$. 
    inline Factory<Zeros, utopia::DYNAMIC> zeros(const Size &size)
    {
        return Factory<Zeros, utopia::DYNAMIC>(size);
    }

    ///nnzXRowOrCol depends if your using a row-major or col-major sparse storage
    template<typename T>
    inline Factory<NNZ<T>, 2> sparse(const Size::SizeType rows, const Size::SizeType cols, T nnzXRowOrCol)
    {
        return Factory<NNZ<T>, 2>(Size({rows, cols}), NNZ<T>(nnzXRowOrCol));
    }


    ///  Returns global matrix \f$ value * 1^{rows \times cols}  \f$.
    template<typename T>
    inline Factory<Values<T>, 2> values(const Size::SizeType rows, const Size::SizeType cols, T value)
    {
        return Factory<Values<T>, 2>(Size({rows, cols}), Values<T>(value));
    }

    ///  Returns global vector \f$ value *1I^{rows \times 1}  \f$.
    template<typename T>
    inline Factory<Values<T>, 1> values(const Size::SizeType rows, T value)
    {
        return Factory<Values<T>, 1>(Size({rows, 1}), Values<T>(value));
    }
     /** @}*/


    /** \addtogroup Local
     * @brief Manipulation with objects on local adress space.
     * @ingroup factory
     *  @{
     */

    /// Returns local identity matrix  \f$ I^{row \times cols}  \f$ i.e. each processors owns local identity matrix. 
    inline Factory<LocalIdentity, 2> local_identity(const Size::SizeType rows, const Size::SizeType cols)
    {
        return Factory<LocalIdentity, 2>(Size({rows, cols}));
    }

    /// Returns local identity matrix  \f$ I^{size \times size}  \f$ i.e. each processors owns local identity matrix. 
    inline Factory<LocalIdentity, 2> local_identity(const Size &size)
    {
        return Factory<LocalIdentity, 2>(size);
    }


    ///  Returns local zero vector \f$ 0^{n \times 1}  \f$. 
    inline Factory<LocalZeros, 1> local_zeros(const Size::SizeType n)
    {
        return Factory<LocalZeros, 1>(Size({n}));
    }

    ///  Returns local zero matrices \f$ 0^{size \times size}  \f$ i.e. each processors owns local zero matrix. 
    inline Factory<LocalZeros, utopia::DYNAMIC> local_zeros(const Size &size)
    {
        return Factory<LocalZeros, utopia::DYNAMIC>(size);
    }

    ///  Returns global matricx \f$ ?^{size_0 \times size_1}  \f$.
    inline Factory<Resize, utopia::DYNAMIC> dense(const Size &size)
    {
        return Factory<Resize, utopia::DYNAMIC>(size);
    }

    ///  Returns global matricx \f$ ?^{size \times size}  \f$.
    inline Factory<Resize, 2> dense(const Size::SizeType rows, const Size::SizeType cols)
    {
        return Factory<Resize, 2>(Size({rows, cols}));
    }


    ///  Returns local matrices \f$ value * 1^{rows \times cols}  \f$ i.e. each processors owns local matrix. 
    template<typename T>
    inline Factory<LocalValues<T>, 2> local_values(const Size::SizeType rows, const Size::SizeType cols, T value)
    {
        return Factory<LocalValues<T>, 2>(Size({rows, cols}), LocalValues<T>(value));
    }

    ///  Returns local vector \f$ value * 1^{rows \times 1}  \f$ i.e. each processors owns local vector. 
    template<typename T>
    inline Factory<LocalValues<T>, 1> local_values(const Size::SizeType rows, T value)
    {
        return Factory<LocalValues<T>, 1>(Size({rows, 1}), LocalValues<T>(value));
    }



    ///nnzXRowOrCol depends if your using a row-major or col-major sparse storage
    template<typename T>
    inline Factory<LocalNNZ<T>, 2> local_sparse(const Size::SizeType rows, const Size::SizeType cols, T nnzXRowOrCol)
    {
        return Factory<LocalNNZ<T>, 2>(Size({rows, cols}), LocalNNZ<T>(nnzXRowOrCol));
    }

    template<typename T>
    inline Factory<LocalRowNNZ<T>, 2> local_row_sparse(const Size::SizeType local_rows, const Size::SizeType global_cols, T nnzXRowOrCol)
    {
        return Factory<LocalRowNNZ<T>, 2>(Size({local_rows, global_cols}), LocalRowNNZ<T>(nnzXRowOrCol));
    }

     /** @}*/

}

#endif //UTOPIA_UTOPIA_FACTORY_HPP
