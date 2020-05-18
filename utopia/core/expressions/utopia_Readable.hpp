//
// Created by Patrick Zulian on 22/05/15.
//

#ifndef UTOPIA_UTOPIA_READABLE_HPP
#define UTOPIA_UTOPIA_READABLE_HPP

#include "utopia_Config.hpp"
#include "utopia_Enums.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    /*!
     * @brief
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template <class Implementation, class Derived, int Order>
    class Readable;

    /*!
     * @brief Specialization of Readable for 2nd order tensors.
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template <class Implementation, class Derived>
    class Readable<Implementation, Derived, 2> {
    public:
        using Scalar = typename Traits<Implementation>::Scalar;

        /**
         * @ingroup     element_acess
         * @brief       Reads value of the element from matrix defined by index (row, column).
         * @warning     Please do not forget to use this function inside of read lock. \n
         *
         * @param[in]  indices  The set of indices.
         * @param[in]  value    The value.
         */
        inline Scalar get(const int row, const int col) const {
            assert_enabled(is_read_locked());
            assert(row < size(derived()).get(0));
            assert(col < size(derived()).get(1));

            return Backend<Scalar, Traits<Implementation>::Backend>::get(derived().implementation(), row, col);
            // return derived().implementation().get(row, col);
        }

    private:
        CONST_DERIVED_CRT(Derived);
    };

    /*!
     * @brief Specialization of Readable for 1st order tensors
     * It provides the interfaces for reading the entries of
     * Implementation.
     * @tparam Implementation the backend type
     * @tparam Derived the derived expression
     */
    template <class Implementation, class Derived>
    class Readable<Implementation, Derived, 1> {
    public:
        using Scalar = typename Traits<Implementation>::Scalar;

        /**
         * @ingroup     element_acess
         * @brief       Gets value of the element, which index matches with requested one.
         * @warning     Please do not forget to use this function inside of read lock. \n
         *
         * @param[in]  indices  The set of indices.
         * @param[in]  value    The value.
         */
        inline Scalar get(const int index) const {
            assert_enabled(is_read_locked());
            assert(index < size(derived()).get(0));
            return Backend<Scalar, Traits<Implementation>::Backend>::get(derived().implementation(), index);
            // return derived().implementation()[index];
        }

        template <typename I, typename T>
        inline void get(const std::vector<I> &index, std::vector<T> &values) const {
            Backend<Scalar, Traits<Implementation>::Backend>::get(derived().implementation(), index, values);
        }

    private:
        CONST_DERIVED_CRT(Derived);
    };

    template <class Tensor>
    class Read {
    public:
        using Scalar = typename Traits<Tensor>::Scalar;

        /**
         * @ingroup     lock
         * @brief       Reading lock providing memory access to the object.
         * @param      tensor  The tensor (the const qualifier might be removed internally).
         */
        Read(const Tensor &tensor) : tensor_(tensor) { const_cast<Tensor &>(tensor_).read_lock(); }

        ~Read() { const_cast<Tensor &>(tensor_).read_unlock(); }

    private:
        const Tensor &tensor_;
    };

    template <class T, int Order>
    class Read<std::vector<Tensor<T, Order> > > {
    public:
        using Tensor = utopia::Tensor<T, Order>;
        using Tensors = std::vector<Tensor>;
        using Scalar = typename Traits<T>::Scalar;

        Read(const Tensors &tensors) : tensors_(tensors) {
            for (auto &t : tensors_) {
                const_cast<Tensor &>(t).read_lock();
            }
        }

        ~Read() {
            for (auto &t : tensors_) {
                const_cast<Tensor &>(t).read_unlock();
            }
        }

        const Tensors &tensors_;
    };

    template <class Tensor>
    class ReadAndWrite {
    public:
        using Scalar = typename Traits<Tensor>::Scalar;

        /**
         * @ingroup     lock
         * @brief       Lock providing both: read and write memory access to the object.
         * @param      tensor  The tensor.
         */
        ReadAndWrite(Tensor &tensor, WriteMode mode = utopia::LOCAL) : tensor_(tensor), mode_(mode) {
            tensor_.read_and_write_lock(mode_);
        }

        ~ReadAndWrite() { tensor_.read_and_write_unlock(mode_); }

    private:
        Tensor &tensor_;
        WriteMode mode_;
    };

}  // namespace utopia
#endif  // UTOPIA_UTOPIA_READABLE_HPP
