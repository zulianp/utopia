//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_WRITABLE_HPP
#define UTOPIA_UTOPIA_WRITABLE_HPP

#include "utopia_Base.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Config.hpp"

#ifdef WITH_PETSC
#include "utopia_PETScTraits.hpp"
#endif //WITH_PETSC

namespace utopia {
    template<class Implementation, class Derived, int Order>
    class Writeable;

    template<class Implementation, class Derived>
    class Writeable<Implementation, Derived, 2> {
    public:
        typedef typename Traits<Implementation>::Scalar Scalar;
        typedef typename Traits<Implementation>::SizeType SizeType;

        /** @defgroup read_write Read/Write
         *  @brief Reading and writing from/to and object
         *  \sa \ref edsl
         */

        /** @defgroup lock Locks
         * @ingroup read_write
         *  @brief  Provides synchronized access to a shared object. 
         */


        /**
         *@ingroup element_acess
         * 
         * @brief       Sets element of matrix on position (row, column) to given value. 
         * @warning     Please do not forget to use this function inside of write lock.  \n
         *              Petsc does not allow to mix add and set.
         * 
         * @param[in]  row      The row index.
         * @param[in]  col      The column index. 
         * @param[in]  value    The value.
         */
        inline void set(const SizeType row, const SizeType col, const Scalar value)
        {
            assert_enabled(is_write_locked());
            assert(row < size(derived()).get(0));
            assert(col < size(derived()).get(1));
            
            Backend<Scalar, Traits<Implementation>::Backend >::Instance().set(derived().implementation(), row, col, value);
        }

        /**
         * @ingroup     element_acess
         * @brief       Adds value to already existing value of the element on position (row, column). 
         * @warning     Please do not forget to use this function inside of write lock. \n
         *              Petsc does not allow to mix add and set.
         * 
         * @param[in]  row      The row index.
         * @param[in]  col      The column index. 
         * @param[in]  value    The value.
         */
        inline void add(const SizeType row, const SizeType col, const Scalar value)
        {
            assert_enabled(is_write_locked());
            assert(row < size(derived()).get(0));
            assert(col < size(derived()).get(1));

            Backend<Scalar, Traits<Implementation>::Backend >::Instance().add(derived().implementation(), row, col, value);
        }


        /**
         * @ingroup     element_acess
         * @brief       Sets prescribed value to all elements, which index matches with one provided by index sets. 
         * @warning     Please do not forget to use this function inside of write lock. \n
         * 
         * @param[in]  row      The set of row indices.
         * @param[in]  col      The set of column indices. 
         * @param[in]  value    The value.
         */
        template<typename Ordinal>
        inline void set(const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend >::Instance().set(derived().implementation(), rows, columns, values);
        }

        template<typename RowT, typename ColT, typename ScalarT>
        inline void set(std::initializer_list<RowT> rows, std::initializer_list<ColT> cols, std::initializer_list<ScalarT> values)
        {
            assert_enabled(is_write_locked());

            using std::copy;
            std::vector<SizeType> vrows(rows.size()), vcols(cols.size());
            std::vector<Scalar> vvalues(values.size());
            copy(rows.begin(), rows.end(), vrows.begin());
            copy(cols.begin(), cols.end(), vcols.begin());
            copy(values.begin(), values.end(), vvalues.begin());

            set(vrows, vcols, vvalues);
        }

#ifdef ENABLE_LOCK_CHECK
        Writeable()
        : lock_active_(false)
        { }

        inline bool is_write_locked() const
        {
            return lock_active_;
        }

        inline void write_lock()
        {
            lock_active_ = true;
        }

        inline void write_unlock()
        {
            lock_active_ = false;
        }
#endif //ENABLE_LOCK_CHECK   

    private:
        DERIVED_CRT(Derived);
#ifdef ENABLE_LOCK_CHECK
        bool lock_active_;
#endif //ENABLE_LOCK_CHECK  
    };

    template<class Implementation, class Derived>
    class Writeable<Implementation, Derived, 1> {
    public:
        typedef typename Traits<Implementation>::Scalar Scalar;

        /**
         * @ingroup     element_acess
         * @brief       Sets element of vector on position (i) to given value. 
         * @warning     Please do not forget to use this function inside of write lock.  \n
         *              Petsc does not allow to mix add and set.
         * @param[in]  index  The index.
         * @param[in]  value  The value.
         */
        inline void set(const SizeType index, const Scalar value)
        {
            assert_enabled(is_write_locked());
            assert(index < size(derived()).get(0));
            Backend<Scalar, Traits<Implementation>::Backend >::Instance().set(derived().implementation(), index, value);
        }

        /**
         * @ingroup     element_acess
         * @brief       Adds prescribed value to the value of the element on the index i.
         * @warning     Please do not forget to use this function inside of write lock. \n
         *              Petsc does not allow to mix add and set.
         * @param[in]  index  The index.
         * @param[in]  value  The value.
         */
        inline void add(const SizeType index, const Scalar value)
        {
            assert_enabled(is_write_locked());
            assert(index < size(derived()).get(0));
            Backend<Scalar, Traits<Implementation>::Backend >::Instance().add(derived().implementation(), index, value);
        }


        /**
         * @ingroup     element_acess
         * @brief       Sets prescribed value to all elements, which index matches with one provided by index set. 
         * @warning     Please do not forget to use this function inside of write lock. \n
         * 
         * @param[in]  indices  The set of indices. 
         * @param[in]  value    The value. 
         */
        template<typename Ordinal>
        inline void set(const std::vector<Ordinal> &indices, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend >::Instance().set(derived().implementation(), indices, values);
        }

#ifdef ENABLE_LOCK_CHECK
        Writeable()
        : lock_active_(false)
        {
        }

        inline bool is_write_locked() const
        {
            return lock_active_;
        }

        inline void write_lock()
        {
            lock_active_ = true;
        }

        inline void write_unlock()
        {
            lock_active_ = false;
        }
#endif //ENABLE_LOCK_CHECK   

    private:
        DERIVED_CRT(Derived);
#ifdef ENABLE_LOCK_CHECK
        bool lock_active_;
#endif //ENABLE_LOCK_CHECK  
    };


    template<class Tensor>
    class Write {
    public:
        Tensor &_tensor;
        typedef typename Traits<Tensor>::Scalar Scalar;

        /**
         * @ingroup     lock
         * @brief       Write lock providing memory access to object. \n
         * @param       tensor  The tensor.
         */
        Write(Tensor &tensor)
        : _tensor(tensor)
        {
#ifdef ENABLE_LOCK_CHECK
            _tensor.write_lock();
#endif //ENABLE_LOCK_CHECK            
            Backend<Scalar, Traits<Tensor>::Backend >::Instance().writeLock(_tensor.implementation());
        }

        ~Write()
        {
#ifdef ENABLE_LOCK_CHECK
            _tensor.write_unlock();
#endif //ENABLE_LOCK_CHECK   
            Backend<Scalar, Traits<Tensor>::Backend >::Instance().writeUnlock(_tensor.implementation());
        }
    };
}

#endif //UTOPIA_UTOPIA_WRITABLE_HPP
