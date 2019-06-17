//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_WRITABLE_HPP
#define UTOPIA_UTOPIA_WRITABLE_HPP

#include "utopia_Base.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Config.hpp"

#include <vector>
#include <initializer_list>

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
            assert(row < (SizeType)size(derived()).get(0));
            assert(col < (SizeType)size(derived()).get(1));

            Backend<Scalar, Traits<Implementation>::Backend>::set(derived().implementation(), row, col, value);
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
            assert(row < (SizeType)size(derived()).get(0));
            assert(col < (SizeType)size(derived()).get(1));

            Backend<Scalar, Traits<Implementation>::Backend>::add(derived().implementation(), row, col, value);
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
        // template<typename Ordinal>
        // inline void set(const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns, const std::vector<Scalar> &values)
        // {
        //     assert_enabled(is_write_locked());
        //     Backend<Scalar, Traits<Implementation>::Backend>::set(derived().implementation(), rows, columns, values);
        // }

        template<typename Ordinal>
        inline void add_matrix(const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend>::add_matrix(derived().implementation(), rows, columns, values);
        }

        inline void set_matrix(const std::vector<SizeType> &rows, const std::vector<SizeType> &columns, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend>::set_matrix(derived().implementation(), rows, columns, values);
        }

        template<typename Ordinal>
        inline void set_matrix(const std::vector<Ordinal> &rows, const std::vector<Ordinal> &columns, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend>::set_matrix(derived().implementation(), rows, columns, values);
        }

        // template<class Rows, class Cols, class Values>
        // inline void set_matrix(const Rows &rows, Cols &cols, const Values &values)
        // {
        //     assert_enabled(is_write_locked());
        //     // Backend<Scalar, Traits<Implementation>::Backend>::set_matrix(derived().implementation(), rows, columns, values);

        //     auto n_rows = rows.size();
        //     auto n_cols = cols.size();

        //     for(std::size_t i = 0; i < n_rows; ++i) {
        //         const auto i_offset = i * n_cols;
        //         for(std::size_t j = 0; j < n_cols; ++j) {
        //             this->set(rows[i], cols[j], values[i_offset + j]);
        //         }
        //     }
        // }

        
        template<class Rows, class Cols, class Values>
        inline void set_matrix(const Rows &rows, Cols &cols, const Values &values)
        {
            assert_enabled(is_write_locked());

            auto n_rows = rows.size();
            auto n_cols = cols.size();

            std::vector<SizeType> ir(1), ic(1);
            std::vector<Scalar> v(1);

            for(std::size_t i = 0; i < n_rows; ++i) {
                const auto i_offset = i * n_cols;
                for(std::size_t j = 0; j < n_cols; ++j) {
                    ir[0] = rows[i];
                    ic[0] = cols[i];
                    v[0]  = values[i_offset + j];
                   
                    Backend<Scalar, Traits<Implementation>::Backend>::set_matrix(
                        derived().implementation(),
                        ir,
                        ic,
                        v
                    );
                }
            }
        }

        template<typename RowT, typename ColT, typename ScalarT>
        inline void set_matrix(std::initializer_list<RowT> rows, std::initializer_list<ColT> cols, std::initializer_list<ScalarT> values)
        {
            assert_enabled(is_write_locked());

            using std::copy;
            std::vector<SizeType> vrows(rows.size()), vcols(cols.size());
            std::vector<ScalarT> vvalues(values.size());
            copy(rows.begin(), rows.end(), vrows.begin());
            copy(cols.begin(), cols.end(), vcols.begin());
            copy(values.begin(), values.end(), vvalues.begin());

            set_matrix(vrows, vcols, vvalues);
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
            assert(index < (SizeType)size(derived()).get(0));
            Backend<Scalar, Traits<Implementation>::Backend>::set(derived().implementation(), index, value);
        }

        inline void set(const Scalar value)
        {
            Backend<Scalar, Traits<Implementation>::Backend>::set(derived().implementation(), value);
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
            assert(index < (SizeType)size(derived()).get(0));
            Backend<Scalar, Traits<Implementation>::Backend>::add(derived().implementation(), index, value);
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
            Backend<Scalar, Traits<Implementation>::Backend>::set(derived().implementation(), indices, values);
        }

        template<typename Ordinal>
        inline void add(const std::vector<Ordinal> &indices, const std::vector<Scalar> &values)
        {
            assert_enabled(is_write_locked());
            Backend<Scalar, Traits<Implementation>::Backend>::add(derived().implementation(), indices, values);
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

    enum WriteMode {
        AUTO  = 0, //unsafe for the moment depends on backend implementation
        LOCAL = 1,
        GLOBAL_INSERT = 2,
        GLOBAL_ADD    = 3
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
        Write(Tensor &tensor, WriteMode mode = utopia::AUTO)
        : _tensor(tensor), mode_(mode)
        {
#ifdef ENABLE_LOCK_CHECK
            _tensor.write_lock();
#endif //ENABLE_LOCK_CHECK
            Backend<Scalar, Traits<Tensor>::Backend>::write_lock(_tensor.implementation(), mode_);
        }

        ~Write()
        {
#ifdef ENABLE_LOCK_CHECK
            _tensor.write_unlock();
#endif //ENABLE_LOCK_CHECK
            Backend<Scalar, Traits<Tensor>::Backend>::write_unlock(_tensor.implementation(), mode_);
        }


        WriteMode mode_;
    };

    template<class T, int Order>
    class Write< std::vector<Wrapper<T, Order> > > {
    public:
        using Tensors = std::vector<Wrapper<T, Order> >;
        using Scalar  = typename Traits<T>::Scalar;
        

        Write(Tensors &tensors, WriteMode mode = utopia::AUTO)
        : tensors_(tensors), mode_(mode)
        {
            for(auto &t : tensors_) {
                Backend<Scalar, Traits<T>::Backend>::write_lock(t.implementation(), mode_);
            }
        }

        ~Write()
        {
            for(auto &t : tensors_) {
                Backend<Scalar, Traits<T>::Backend>::write_unlock(t.implementation(), mode_);
            }
        }

        Tensors &tensors_;
        WriteMode mode_;
    };

}

#endif //UTOPIA_UTOPIA_WRITABLE_HPP
