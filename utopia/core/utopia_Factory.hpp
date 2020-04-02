#ifndef UTOPIA_UTOPIA_FACTORY_HPP
#define UTOPIA_UTOPIA_FACTORY_HPP

#include "utopia_Expression.hpp"
#include "utopia_Size.hpp"
#include "utopia_Optional.hpp"

namespace utopia {

    template<class Type>
    class FactoryTraits {
    public:
        typedef double Scalar;

        static constexpr const char * get_class()
        {
            return "TODO";
        }

        enum {
            FILL_TYPE = FillType::DELEGATE
        };
    };

    class Identity {};
    class DenseIdentity {};

    template<>
    class FactoryTraits<Identity> {
    public:
        typedef double Scalar;

        static constexpr const char * get_class()
        {
            return "Identity";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    template<>
    class FactoryTraits<DenseIdentity> {
    public:
        typedef double Scalar;

        static constexpr const char * get_class()
        {
            return "DenseIdentity";
        }

        enum {
            FILL_TYPE = FillType::DENSE
        };
    };

    class LocalIdentity {};
    class LocalDenseIdentity {};
    class Zeros {};
    class LocalZeros {};


    template<>
    class FactoryTraits<LocalIdentity> {
    public:
        typedef double Scalar;

        static constexpr const char * get_class()
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

        static constexpr const char * get_class()
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

        static constexpr const char * get_class()
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

        static constexpr const char * get_class()
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

        static constexpr const char * get_class()
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

        static constexpr const char * get_class()
        {
            return "NNZ";
        }

        enum {
            FILL_TYPE = FillType::SPARSE
        };
    };

    template<typename _SizeType, typename _IntType, typename _Scalar>
    class CRS {
    public:
        CRS(_SizeType rowPtr, _IntType cols, _Scalar values) : _rowPtr(rowPtr), _cols(cols), _values(values) {}

        inline const _SizeType& rowPtr() const { return _rowPtr; }
        inline const _IntType& cols() const { return _cols; }
        inline const _Scalar& values() const { return _values; }

    private:
        _SizeType _rowPtr;
        _IntType _cols;
        _Scalar _values;
    };

    template<typename SizeType>
    class NNZXRow {
    public:
        NNZXRow(const std::vector<SizeType> &d_nnz, const std::vector<SizeType> &o_nnz)
        : d_nnz(d_nnz), o_nnz(o_nnz)
        {}

        const std::vector<SizeType> &d_nnz, &o_nnz;
    };

    template<typename SizeType>
    class FactoryTraits< NNZXRow<SizeType> > {
    public:
        typedef double Scalar;

        static constexpr const char * get_class()
        {
            return "NNZXRow";
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

        static constexpr const char * get_class()
        {
            return "LocalNNZ";
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

        static constexpr const char * get_class()
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

        inline std::string get_class() const override
        {
            return "Factory(" + std::string(FactoryTraits<Type>::get_class()) + ")";
        }

        virtual ~Factory() {}


    private:
        Size _size;
        Type _type;
    };

    template<class SType, int Order, class Right>
    class MostDescriptive<Factory<SType, Order>, Right> {
    public:
        typedef Right Type;
    };

    template<class Left, class SType, int Order>
    class MostDescriptive<Left, Factory<SType, Order>> {
    public:
        typedef Left Type;
    };


    template<class SType, int Order, class Right>
    class MostDescriptive<Factory<SType, Order>, Number<Right> > {
    public:
        typedef utopia::Factory<SType, Order> Type;
    };

    template<class Left, class SType, int Order>
    class MostDescriptive<Number<Left>, Factory<SType, Order>> {
    public:
        typedef utopia::Factory<SType, Order> Type;
    };

    template<class Factory, class Options>
    class Build : public Expression< Build<Factory, Options> > {
    public:
        static const int Order = Factory::Order;

        enum {
            StoreAs = UTOPIA_BY_VALUE
        };

        typedef typename Factory::Scalar Scalar;

        Build(const Factory &factory, const Options &opts)
        : factory_(factory), opts_(opts)
        {}

        inline std::string get_class() const override
        {
            return factory_.get_class();
        }

        virtual ~Build() {}

        const Options & opts() const
        {
            return opts_;
        }

        const Factory &factory() const
        {
            return factory_;
        }

    private:
        Factory factory_;
        Options opts_;
    };



    template<class Type, int Order_>
    class SymbolicTensor : public Expression< SymbolicTensor<Type, Order_> > {
    public:
        static const int Order = Order_;

        enum {
            StoreAs = UTOPIA_BY_VALUE
        };

        typedef typename FactoryTraits<Type>::Scalar Scalar;

        static inline Type type()
        {
            return Type();
        }

        inline std::string get_class() const override
        {
            return "SymbolicTensor(" + std::string(FactoryTraits<Type>::get_class()) + ")";
        }
    };

    template<class SType, int Order, class Right>
    class MostDescriptive<SymbolicTensor<SType, Order>, Right > {
    public:
        typedef Right Type;
    };

    template<class Left, class SType, int Order>
    class MostDescriptive<Left, SymbolicTensor<SType, Order>> {
    public:
        typedef Left Type;
    };


    template<class SType, int Order, class Right>
    class MostDescriptive<SymbolicTensor<SType, Order>, Number<Right> > {
    public:
        typedef utopia::SymbolicTensor<SType, Order> Type;
    };

    template<class Left, class SType, int Order>
    class MostDescriptive<Number<Left>, SymbolicTensor<SType, Order>> {
    public:
        typedef utopia::SymbolicTensor<SType, Order> Type;
    };

    template<class SType,
             int Order,
             class Right,
             class Default,
             int SparsityLeft,
             int SparsityRight>
    class ChooseType<SymbolicTensor<SType, Order>, Right, Default, SparsityLeft, SparsityRight> {
    public:
        typedef Right Type;
    };

    template<class SType,
             int Order,
             class Left,
             class Default,
             int SparsityLeft,
             int SparsityRight>
    class ChooseType<Left, SymbolicTensor<SType, Order>, Default, SparsityLeft, SparsityRight> {
    public:
        typedef Left Type;
    };


    template<class SType,
             int Order,
             class Right,
             class Default,
             int SparsityLeft,
             int SparsityRight>
    class ChooseType<Factory<SType, Order>, Right, Default, SparsityLeft, SparsityRight> {
    public:
        typedef Right Type;
    };

    template<class SType,
             int Order,
             class Left,
             class Default,
             int SparsityLeft,
             int SparsityRight>
    class ChooseType<Left, Factory<SType, Order>, Default, SparsityLeft, SparsityRight> {
    public:
        typedef Left Type;
    };

    template<class Index>
    class Ghosts : public Expression< Ghosts<Index> > {
    public:

        static const int Order = 1;

        Ghosts(const Size::SizeType &local_size, const Size::SizeType &global_size, const Index &index)
        : local_size_(local_size), global_size_(global_size), index_(index)
        {}

        // Ghosts(const Size::SizeType &local_size, const Size::SizeType &global_size, Index &&index)
        // : local_size_(local_size), global_size_(global_size), index_(std::move(index))
        // {}

        const Index &index() const
        {
            return index_;
        }

        const Size::SizeType &local_size() const
        {
            return local_size_;
        }

        const Size::SizeType &global_size() const
        {
            return global_size_;
        }

    private:
        Size::SizeType local_size_;
        Size::SizeType global_size_;
        Index index_;
    };

    template<class Index>
    class Traits< Ghosts<Index> > {
    public:
        typedef double Scalar;

        enum {
            FILL_TYPE = FillType::DENSE
        };
    };

    template<class Type, int Order>
    class Traits< SymbolicTensor<Type, Order> > {
    public:
        typedef typename utopia::FactoryTraits<Type>::Scalar Scalar;

        enum {
            FILL_TYPE = FactoryTraits<Type>::FILL_TYPE
        };
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

    ///FIXME make proper version
    template<typename _SizeType, typename _IntType, typename _Scalar>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Factory<CRS<_SizeType, _IntType, _Scalar>, 2> crs(const Size::SizeType rows, const Size::SizeType cols, _SizeType &rowPtr, _IntType &crs_columns, _Scalar &values)
    {
        return Factory<CRS<_SizeType, _IntType, _Scalar>, 2>(Size({rows, cols}), CRS<_SizeType, _IntType, _Scalar>(rowPtr, crs_columns, values));
    }


    //FIXME make proper version
    template<class Index>
    UTOPIA_DEPRECATED_MSG("Use the object variant with 1st argument as the communicator.")
    inline Ghosts<Index> ghosted(
        const Size::SizeType &local_size,
        const Size::SizeType &global_size,
        Index &&index)
    {
        return Ghosts<Index>(local_size, global_size, std::forward<Index>(index));
    }
}


#endif //UTOPIA_UTOPIA_FACTORY_HPP
