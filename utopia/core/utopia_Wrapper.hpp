// @file
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_MATRIXEXPR_HPP
#define utopia_utopia_MATRIXEXPR_HPP


#include "utopia_Expression.hpp"
#include "utopia_Evaluator.hpp"
#include "utopia_Assign.hpp"
#include "utopia_Traits.hpp"
#include "utopia_InPlace.hpp"
#include "utopia_Mutable.hpp"
#include "utopia_Readable.hpp"
#include "utopia_Writable.hpp"
#include "utopia_Ranged.hpp"
#include "utopia_Select.hpp"


#include <iostream>
#include <type_traits>


#define TPL_REMOVE_CONST(T) typename std::remove_const<T>::type


namespace utopia {
    template<int Order>
    struct Order2String {
        static std::string Value() {
            return "Undefined";
        }
    };

    template<>
    struct Order2String<1> {
        static std::string Value() {
            return "Vector";
        }
    };

    template<>
    struct Order2String<2> {
        static std::string Value() {
            return "Matrix";
        }
    };


    /*!
     * @class Wrapper
     * @brief The wrapper of any tensor representation.
     */
    template<class _Implementation, int _Order>
    class Wrapper : public Expression<Wrapper<_Implementation, _Order> >,
                    public Mutable<_Implementation, Wrapper<_Implementation, _Order> >,
                    public Readable<_Implementation, Wrapper<_Implementation, _Order>, _Order>,
                    public Writeable<_Implementation, Wrapper<_Implementation, _Order>, _Order>,
                    public Structured< Wrapper<_Implementation, _Order> >,
                    public Ranged< Wrapper<_Implementation, _Order>, _Order>,
                    public Selectable< _Implementation, Wrapper<_Implementation, _Order>, _Order >
            {
    public:
        typedef _Implementation Implementation;
        typedef typename Traits<Implementation>::Scalar Scalar;
        typedef typename Traits<Implementation>::SizeType SizeType;

        enum {
            Backend = Traits<Implementation>::Backend
        };

        enum {
            FILL_TYPE = Traits<Implementation>::FILL_TYPE
        };

        static const int Order = _Order;

        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        virtual ~Wrapper() {}

        template<class Derived>
        Wrapper(const Expression<Derived> &expr) {
            evaluator().eval(Construct<Wrapper, Derived>(*this, expr.derived()));
        }

        template<class Derived>
        Wrapper &operator=(const Expression<Derived> &expr) {
            evaluator().eval(Assign<Wrapper, Derived>(*this, expr.derived()));
            return *this;
        }

        Wrapper & operator=(Wrapper &&other) {
             utopia::Backend<Scalar, Traits<Implementation>::Backend>::Instance().assign(_impl, std::move(other._impl));
            return *this;
        }

        Wrapper & operator=(const Wrapper &other) {
            utopia::Backend<Scalar, Traits<Implementation>::Backend>::Instance().assign(_impl, other._impl);
            return *this;
        }


        Wrapper(const Wrapper &expr) {
            // evaluator().eval(Construct<Wrapper, Wrapper>(*this, expr));
            //FIXME
            _impl = expr._impl;
        }
        
        Wrapper() { }

        Wrapper(const Implementation &impl)
                : _impl(impl) { }

        Wrapper(Implementation &&impl)
                : _impl(std::move(impl)) { }


        const Implementation &implementation() const {
            return _impl;
        }

        Implementation &implementation() {
            return _impl;
        }


        std::string getClass() const {
            return Fill2String<FILL_TYPE>::Value() + Order2String<_Order>::Value();
        }


        Evaluator<Implementation> &evaluator() { return _evaluator; }

    private:
        Implementation _impl;
        Evaluator<Implementation> _evaluator;
    };

    /*!
     * @brief The wrapper of a reference to any tensor representation.
     */
    template<class _Implementation, int _Order>
    class Wrapper<_Implementation &, _Order> :  public Expression<Wrapper<_Implementation &, _Order> >,
                                                public Mutable<_Implementation, Wrapper<_Implementation &, _Order> >,
                                                public Readable<_Implementation, Wrapper<_Implementation &, _Order>, _Order>,
                                                public Writeable<_Implementation, Wrapper<_Implementation &, _Order>, _Order>,
                                                public Structured< Wrapper<_Implementation &, _Order> >,
                                                public Ranged< Wrapper<_Implementation &, _Order>, _Order>,
                                                public Selectable< _Implementation, Wrapper<_Implementation, _Order>, _Order >

    {
    public:
        typedef _Implementation Implementation;
        typedef typename Traits<Implementation>::Scalar Scalar;
        typedef typename Traits<Implementation>::SizeType SizeType;

        enum {
            Backend = Traits<Implementation>::Backend
        };

        static const int Order = _Order;

        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };

        enum {
            FILL_TYPE = Traits<Implementation>::FILL_TYPE
        };


        virtual ~Wrapper() { }


        template<class Derived>
        Wrapper &operator=(const Expression<Derived> &expr) {
            evaluator().eval(Assign<Wrapper, Derived>(*this, expr.derived()));
            return *this;
        }

        Wrapper(Implementation &impl)
                : _impl(impl) {
            static_assert(!std::is_const<Implementation>::value, "Wrong specialization instantiated");
        }


        const Implementation &implementation() const {
            return _impl;
        }

        Implementation &implementation() {
            return _impl;
        }


       std::string getClass() const {
           return Fill2String<FILL_TYPE>::Value() + Order2String<_Order>::Value();
       }



        Evaluator<Implementation> &evaluator() { return _evaluator; }

    private:
        Implementation &_impl;
        Evaluator<Implementation> _evaluator;
    };


    /*!
     * @brief The wrapper of a const reference to any tensor representation.
     */
    template<class _Implementation, int _Order>
    class Wrapper<const _Implementation &, _Order> : public Expression< Wrapper<const _Implementation &, _Order> >,
                                                     public Readable<_Implementation, Wrapper<const _Implementation &, _Order>, _Order>,
                                                     public Structured< Wrapper<const _Implementation &, _Order> >,
                                                     public Selectable< _Implementation, Wrapper<_Implementation, _Order>, _Order > {
    public:
        typedef _Implementation Implementation;
        typedef typename Traits<Implementation>::Scalar Scalar;
        typedef typename Traits<Implementation>::SizeType SizeType;


        virtual ~Wrapper() { }

        const Implementation &implementation() const {
            return _impl;
        }

        Wrapper(const Implementation &impl)
                : _impl(impl) { }


        static const int Order = _Order;

        enum {
            FILL_TYPE = Traits<Implementation>::FILL_TYPE
        };

        enum {
            StoreAs = UTOPIA_BY_REFERENCE
        };


        std::string getClass() const {
            return Fill2String<FILL_TYPE>::Value() + Order2String<_Order>::Value();
        }
    private:

        const Implementation &_impl;
    };

    template<class Tensor, int Order>
    inline auto backend(const Wrapper<Tensor, Order> &) -> decltype(Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance()) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance();
    }

    template<class Tensor, int Order>
    inline const BackendInfo &backend_info(const Wrapper<Tensor, Order> &t) {
        return backend(t).info();
    }


    template<class Tensor, int Order>
    inline long unique_id(const Wrapper<Tensor, Order> &t)
    {
        //FIXME
        return (long) &t;
    } 


    /** \addtogroup ranges
     * @ingroup read_write
     * @brief Ranges provide the necessary access information for distributed memory data structures such as vectors
     * and matrices.
     *  @{
     */

    /*!
     * @fn Range range(const Wrapper<Tensor, 1> &v) 
     * @brief  ranges allow to coordinate the access of parts of the tensor based on memory location.
     * The location is defined by how the vector is paritioned among processes.
     * @tparam Tensor the backend type of the 1st order tensor
     * @return the range of the input vector v.
     */
    template<class Tensor>
    inline Range range(const Wrapper<Tensor, 1> &v) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().range(v.implementation());
    }

    /*!
     * @fn Range row_range(const Wrapper<Tensor, 2> &v)
     * @brief Row ranges allow to coordinate the access of rows of the tensor based on memory location.
     * The location is defined by how the matrix is paritioned among processes.
     * @tparam Tensor the backend type of the 2nd order tensor
     * @return the row range of the input matrix v.
     */
    template<class Tensor>
    inline Range row_range(const Wrapper<Tensor, 2> &v) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().row_range(
                v.implementation());
    }

    /*!
     * @fn Range col_range(const Wrapper<Tensor, 2> &v)
     * @brief Column ranges allow to coordinate the access of columns of the tensor based on memory location.
     * The location is defined by how the matrix is paritioned among processes.
     * @tparam Tensor the backend type of the 2nd order tensor
     * @return the column range of the input matrix v.
     */
    template<class Tensor>
    inline Range col_range(const Wrapper<Tensor, 2> &v) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().col_range(
                v.implementation());
    }

    /** @}*/



    /**     @defgroup   io Input/Output
     *      @ingroup    base_functions
     */


    /**
     * @ingroup    io
     * @brief      Reads tensor from file.
     *
     * @param[in]  path    The path. 
     * @param[in]  t       Tensor to be read and used. 
     */
    template<class Tensor, int Order>
    inline bool read(const std::string &path, Wrapper<Tensor, Order> &t) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().read(path, t.implementation());
    }

    template<class Tensor, int Order, class... Args>
    inline bool read(const std::string &path, Wrapper<Tensor, Order> &t, Args &&...args) {
    	auto &backend = Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance();
        return backend.read(path, t.implementation(), backend.parse_args(options(args...)));
    }

    /**
     * @ingroup    io
     * @brief      Writes and saves tensor into file.
     *
     * @param[in]  path    The path. 
     * @param[in]  t       Tensor to be saved. 
     */
    template<class Tensor, int Order>
    inline bool write(const std::string &path, const Wrapper<Tensor, Order> &t) {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().write(path, t.implementation());
    }

    /**
     * @defgroup   profiling Profiling/Debugging
     * @ingroup    base_functions
     */


    /**
     * @defgroup   interoperability Interoperability
     * @ingroup    base_functions
     */

    /**
     * @ingroup    profiling
     * @brief      Monitoring of the specific tensor over iterations in solvers. \n
     *             Outputs matlab file, with saved values of tensor over iterates. \n
     *
     * @param[in]  it      The iteration. 
     * @param      t       Tensor to be monitored.
     */
    template<class Tensor, int Order>
    inline void monitor(const long &it, Wrapper<Tensor, Order> &t) 
    {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().monitor(it, t.implementation());
    }

    /**
     * @ingroup    profiling
     * @brief      Monitoring of the specific tensor over iterations in solvers. \n
     *             Outputs matlab file, with saved values of tensor over iterates. \n
     *
     * @param[in]  it                    The iteration number. 
     * @param[in]  name_of_file          The name of file.
     * @param[in]  name_of_instance      The name of tensor to be saved.
     * @param      t                     Tensor to be monitored.
     */
    template<class Tensor, int Order>
    inline void monitor(const long &it, Wrapper<Tensor, Order> &t, const std::string name_of_file, const std::string name_of_instance) 
    {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().monitor(it, t.implementation(), name_of_file, name_of_instance);
    }



    /**
     * @ingroup    profiling
     * @brief      Gets number of global nnz.
     *
     * @return     The global nnz.
     */
    template<class Tensor>
    inline typename Traits<Tensor>::Scalar get_global_nnz(Wrapper<Tensor, 2> &t) 
    {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().get_global_nnz(t.implementation());
    }

    /**
     * @ingroup    profiling
     * @brief      Gets number of local nnz.
     *
     * @return     The local nnz.
     */
    template<class Tensor>
    inline typename Traits<Tensor>::Scalar get_local_nnz(Wrapper<Tensor, 2> &t) 
    {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().get_local_nnz(t.implementation());
    }

    /**
     * @ingroup    interoperability
     * @brief      Converts backend-type tensor into utopia-wrapper specified type of tensor.
     *
     * @param      rawType  The external library tensor (backend supported).
     * @param      t        Utopia tensor(wrapper). 
     */
    template<class RawType, class Tensor, int Order>
    inline void convert(RawType &rawType, Wrapper<Tensor, Order> &t) {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().convert(rawType, t.implementation());
    }
    
    /**
     * @ingroup    interoperability
     * @brief      Converts utopia-wrapper specified type of tensor into original type. 
     *
     * @param      t        Utopia tensor(wrapper). 
     * @param      rawType  The external library tensor (backend supported).
     */
    template<class Tensor, int Order, class RawType>
    inline void convert(Wrapper<Tensor, Order> &t, RawType &rawType) {
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().convert(t.implementation(), rawType);
    }


    template<typename T, int Order>
    void disp(const Wrapper<std::vector<T>, Order> &w, std::ostream &os) {
        disp(w.implementation().begin(), w.implementation().end(), os);
    }

    inline void disp(const double value, std::ostream &os = std::cout)
    {
        os << value << "\n";
    }

    // inline void disp(const double value)
    // {
    //     disp(value, std::cout);
    // }

    template<class VectorT>
    Wrapper<VectorT, 1> vmake() {
        return Wrapper<VectorT, 1>();
    }

    template<class VectorT>
    Wrapper<VectorT, 1> vmake(const VectorT &v) {
        return Wrapper<VectorT, 1>(v);
    }

    template<class VectorT, typename...Args>
    Wrapper<VectorT, 1> vmake(Args &&...args) {
        return Wrapper<VectorT, 1>(VectorT(std::forward<Args>(args)...));
    }

    template<class VectorT, typename T>
    Wrapper<VectorT, 1> vmake(std::initializer_list<T> args) {
        return Wrapper<VectorT, 1>(VectorT(args));
    }

    template<class MatrixT>
    Wrapper<MatrixT, 2> mmake() {
        return Wrapper<MatrixT, 2>();
    }

    template<class MatrixT, typename...Args>
    Wrapper<MatrixT, 2> mmake(Args &&...args) {
        return Wrapper<MatrixT, 2>(MatrixT(std::forward<Args>(args)...));
    }

    template<class VectorT>
    Wrapper<VectorT &, 1> vref(VectorT &v) {
        return v;
    }

    template<class VectorT>
    Wrapper<const VectorT &, 1> vref(const VectorT &v) {
        return v;
    }

    template<typename MatrixT>
    Wrapper<MatrixT &, 2> mref(MatrixT &m) {
        return m;
    }

    template<typename MatrixT>
    Wrapper<const MatrixT &, 2> mref(const MatrixT &m) {
        return m;
    }


    template<int Order, typename Tensor>
    inline Wrapper<Tensor &, Order> wrap(Tensor &t)
    {
        return t;
    }

    template<int Order, typename Tensor>
    inline Wrapper<std::shared_ptr<Tensor>, Order> wrap(const std::shared_ptr<Tensor> &t)
    {
        return t;
    }

    /**
     * @ingroup    io 
     * @brief      Displays any tensor. 
     *
     * @param[in]  w      The tensor tobe displayed. 
     *
     * @tparam     Impl   The tensor implementation.
     * @tparam     Order  The order of the tensor. 
     */
    template<class Impl, int Order>
    void disp(const Wrapper<Impl, Order> &w)
    {
        // disp(w, std::cout);
        return Backend<typename Traits<Impl>::Scalar, Traits<Impl>::Backend>::Instance().disp(w.implementation());
    }

    /**
     * @ingroup    queries
     * @brief      Checks, if Tensor was assembled. 
     *
     * @param[in]  w               The wrapper of Tensor. 
     *
     * @tparam     Implementation  Actual tensor. 
     * @tparam     Order           The order of the tensor. 
     *
     * @return     State of assembly. 
     */
    template<class Implementation, int Order>
    bool empty(const Wrapper<Implementation, Order> &w)
    {
        static_assert(Order >= 1, "Does not work for scalars");
        auto s = size(w); 
        return s.get(0) == INVALID_INDEX;
    }

    /**
     * @ingroup    queries
     * @brief      Checks, if Tensor contains inf/nan. 
     *
     * @param[in]  w       The wrapper of Tensor. 
     *
     * @tparam     Tensor  Actual tensor. 
     * @tparam     Order   The order of the tensor. 
     *
     * @return     1 if there is some nan or inf, 0 otherwise
     */
    template<class Tensor, int Order>
    bool has_nan_or_inf(const Wrapper<Tensor, Order> &w)
    {
        // static_assert(Order == 1, "contains_nan_or_inf:: works just for vectors at the moment");
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().is_nan_or_inf(w.implementation());
    }

    template<class Derived>
    inline constexpr int order(const Expression<Derived> &)
    {
        return Derived::Order;
    }

    template<class Tensor, int Order>
    inline auto raw_type(const Wrapper<Tensor, Order> &w) -> decltype( 
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().raw_type(w.implementation()) 
        ) &
    {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().raw_type(w.implementation());
    }

    template<class Tensor, int Order>
    inline auto raw_type(Wrapper<Tensor, Order> &w) -> decltype( 
        Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().raw_type(w.implementation())
        ) &
    {
        return Backend<typename Traits<Tensor>::Scalar, Traits<Tensor>::Backend>::Instance().raw_type(w.implementation());
    }
}

#endif //utopia_utopia_MATRIXEXPR_HPP
