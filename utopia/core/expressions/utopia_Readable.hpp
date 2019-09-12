//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_READABLE_HPP
#define UTOPIA_UTOPIA_READABLE_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Config.hpp"
#include "utopia_Enums.hpp"

namespace utopia {

    /*!
     * @brief
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template<class Implementation, class Derived, int Order>
    class Readable;

    /*!
     * @brief Specialization of Readable for 2nd order tensors.
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template<class Implementation, class Derived>
    class Readable<Implementation, Derived, 2> {
    public:
        typedef typename Traits<Implementation>::Scalar Scalar;

         /**
         * @ingroup     element_acess
         * @brief       Reads value of the element from matrix defined by index (row, column).
         * @warning     Please do not forget to use this function inside of read lock. \n
         *
         * @param[in]  indices  The set of indices.
         * @param[in]  value    The value.
         */
        inline Scalar get(const int row, const int col) const
        {
            assert_enabled(is_read_locked());
            assert(row < size(derived()).get(0));
            assert(col < size(derived()).get(1));

            return Backend<Scalar, Traits<Implementation>::Backend >::get(derived().implementation(), row, col);
            // return derived().implementation().get(row, col);
        }

#ifdef ENABLE_LOCK_CHECK
        Readable()
        : lock_active_(false)
        { }

        inline bool is_read_locked() const
        {
            return lock_active_;
        }

        inline void read_lock() const
        {
            lock_active_ = true;
        }

        inline void read_unlock() const
        {
            lock_active_ = false;
        }
#endif //NDEBUG

    private:
        CONST_DERIVED_CRT(Derived);

#ifdef ENABLE_LOCK_CHECK
        mutable bool lock_active_;
#endif //NDEBUG

    };

    /*!
     * @brief Specialization of Readable for 1st order tensors
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template<class Implementation, class Derived>
    class Readable<Implementation, Derived, 1> {
    public:
        typedef typename Traits<Implementation>::Scalar Scalar;

         /**
         * @ingroup     element_acess
         * @brief       Gets value of the element, which index matches with requested one.
         * @warning     Please do not forget to use this function inside of read lock. \n
         *
         * @param[in]  indices  The set of indices.
         * @param[in]  value    The value.
         */
        inline Scalar get(const int index) const
        {
            assert_enabled(is_read_locked());
            assert(index < size(derived()).get(0));
            return Backend<Scalar, Traits<Implementation>::Backend >::get(derived().implementation(), index);
            // return derived().implementation()[index];

        }

        template<typename I, typename  T>
        inline void get(const std::vector<I> &index, std::vector<T> &values) const
        {
            Backend<Scalar, Traits<Implementation>::Backend >::get(derived().implementation(), index, values);
        }

// #ifdef ENABLE_LOCK_CHECK
//         Readable()
//         : lock_active_(false)
//         { }

//         inline bool is_read_locked() const
//         {
//             return lock_active_;
//         }

//         inline void read_lock() const
//         {
//             lock_active_ = true;
//         }

//         inline void read_unlock() const
//         {
//             lock_active_ = false;
//         }
// #endif //NDEBUG

    private:
        CONST_DERIVED_CRT(Derived);

// #ifdef ENABLE_LOCK_CHECK
//         mutable bool lock_active_;
// #endif //NDEBUG
    };


    template<class Tensor>
    class Read {
    public:
        typedef typename Traits<Tensor>::Scalar Scalar;


        /**
         * @ingroup     lock
         * @brief       Reading lock providing memory access to the object.
         * @param      tensor  The tensor (the const qualifier might be removed internally).
         */
        Read(const Tensor &tensor)
        : tensor_(tensor)
        {
            const_cast<Tensor &>(tensor_).read_lock();
        }

        ~Read()
        {
            const_cast<Tensor &>(tensor_).read_unlock();
        }

    private:
        const Tensor &tensor_;
    };


    template<class T, int Order>
    class Read< std::vector<Tensor<T, Order> > > {
    public:
        using Tensor = utopia::Tensor<T, Order>;
        using Tensors = std::vector<Tensor>;
        using Scalar  = typename Traits<T>::Scalar;
        

        Read(const Tensors &tensors)
        : tensors_(tensors)
        {
            for(auto &t : tensors_) {
                const_cast<Tensor &>(t).read_lock();
            }
        }

        ~Read()
        {
            for(auto &t : tensors_) {
                const_cast<Tensor &>(t).read_unlock();
            }
        }

        const Tensors &tensors_;
    };

    template<class Tensor>
    class ReadAndWrite {
    public:
        typedef typename Traits<Tensor>::Scalar Scalar;

        /**
         * @ingroup     lock
         * @brief       Lock providing both: read and write memory access to the object.
         * @param      tensor  The tensor.
         */
        ReadAndWrite(Tensor &tensor, WriteMode mode)
        : tensor_(tensor), mode_(mode)
        {
            tensor_.read_and_write_lock(mode_);
        }

        ~ReadAndWrite()
        {
            tensor_.read_and_write_unlock(mode_);
        }

    private:
         Tensor &tensor_;
         WriteMode mode_;
    };

}
#endif //UTOPIA_UTOPIA_READABLE_HPP
