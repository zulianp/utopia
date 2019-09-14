//
// Created by Patrick Zulian on 26/05/15.
//

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

        static constexpr const char * getClass()
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

        static constexpr const char * getClass()
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

        static constexpr const char * getClass()
        {
            return "DenseIdentity";
        }

        enum {
            FILL_TYPE = FillType::DENSE
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

        static constexpr const char * getClass()
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

        static constexpr const char * getClass()
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

        inline std::string getClass() const override
        {
            return "Factory(" + std::string(FactoryTraits<Type>::getClass()) + ")";
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

        inline std::string getClass() const override
        {
            return factory_.getClass();
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

        inline std::string getClass() const override
        {
            return "SymbolicTensor(" + std::string(FactoryTraits<Type>::getClass()) + ")";
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


    /// Returns identity matrix  \f$ I^{row \times cols}  \f$.
    inline Factory<DenseIdentity, 2> dense_identity(const Size::SizeType rows, const Size::SizeType cols)
    {
        return Factory<DenseIdentity, 2>(Size({rows, cols}));
    }
    /// Returns denDensese_identity matrix  \f$ I^{size_0 \times size_1}  \f$.
    inline Factory<DenseIdentity, 2> dense_identity(const Size &size)
    {
        return Factory<DenseIdentity, 2>(size);
    }

    /// Returns identity matrix  \f$ I^{row \times cols}  \f$.
    inline constexpr SymbolicTensor<Identity, 2> identity()
    {
        return SymbolicTensor<Identity, 2>();
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

    ///nnz_x_row_or_col depends if your using a row-major or col-major sparse storage
    template<typename T>
    inline Factory<NNZ<T>, 2> sparse(const Size::SizeType rows, const Size::SizeType cols, T nnz_x_row_or_col)
    {
        return Factory<NNZ<T>, 2>(Size({rows, cols}), NNZ<T>(nnz_x_row_or_col));
    }

    template<typename SizeType>
    inline Factory<NNZXRow<SizeType>, 2> sparse(
        const Size &gs,
        const std::vector<SizeType> &d_nnz,
        const std::vector<SizeType> &o_nnz)
    {
        return Factory<NNZXRow<SizeType>, 2>(gs, NNZXRow<SizeType>(d_nnz, o_nnz));
    }

    template<typename _SizeType, typename _IntType, typename _Scalar>
    inline Factory<CRS<_SizeType, _IntType, _Scalar>, 2> crs(const Size::SizeType rows, const Size::SizeType cols, _SizeType &rowPtr, _IntType &crs_columns, _Scalar &values)
    {
        return Factory<CRS<_SizeType, _IntType, _Scalar>, 2>(Size({rows, cols}), CRS<_SizeType, _IntType, _Scalar>(rowPtr, crs_columns, values));
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



    ///nnz_x_row_or_col depends if your using a row-major or col-major sparse storage
    template<typename T>
    inline Factory<LocalNNZ<T>, 2> local_sparse(const Size::SizeType rows, const Size::SizeType cols, T nnz_x_row_or_col)
    {
        return Factory<LocalNNZ<T>, 2>(Size({rows, cols}), LocalNNZ<T>(nnz_x_row_or_col));
    }

    template<typename T>
    inline Factory<LocalNNZ<T>, 2> local_sparse(const Size &s, T nnz_x_row_or_col)
    {
        return Factory<LocalNNZ<T>, 2>(s, LocalNNZ<T>(nnz_x_row_or_col));
    }

    template<typename T, class... Args>
    inline auto local_sparse(
        const Size::SizeType rows,
        const Size::SizeType cols,
        T nnz_x_row_or_col,
        Args &&... opts) -> Build< Factory<LocalNNZ<T>, 2>, decltype(options(opts...))>
    {
        return Build<Factory<LocalNNZ<T>, 2>, decltype(options(opts...))>(local_sparse(rows, cols, nnz_x_row_or_col), options(opts...));
    }

     /** @}*/

    template<class Index>
    inline Ghosts<Index> ghosted(
        const Size::SizeType &local_size,
        const Size::SizeType &global_size,
        Index &&index)
    {
        return Ghosts<Index>(local_size, global_size, std::forward<Index>(index));
    }
}

#endif //UTOPIA_UTOPIA_FACTORY_HPP
