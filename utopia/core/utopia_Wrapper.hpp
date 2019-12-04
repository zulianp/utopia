#ifndef utopia_utopia_MATRIXEXPR_HPP
#define utopia_utopia_MATRIXEXPR_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Utils.hpp"

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
    struct Order2String<0> {
        static std::string Value() {
            return "Scalar";
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

    // template<class Tensor, int Order>
    // inline long unique_id(const Tensor<Derived, Order> &t)
    // {
    //     //FIXME
    //     return (long) &t;
    // }


    /** \addtogroup ranges
     * @ingroup read_write
     * @brief Ranges provide the necessary access information for distributed memory data structures such as vectors
     * and matrices.
     *  @{
     */

    /*!
     * @fn Range range(const Tensor<Derived, 1> &v)
     * @brief  ranges allow to coordinate the access of parts of the tensor based on memory location.
     * The location is defined by how the vector is paritioned among processes.
     * @tparam Tensor the backend type of the 1st order tensor
     * @return the range of the input vector v.
     */
    template<class Derived>
    inline Range range(const Tensor<Derived, 1> &v) {
        return v.derived().range();
    }

    /*!
     * @fn Range row_range(const Tensor<Derived, 2> &v)
     * @brief Row ranges allow to coordinate the access of rows of the tensor based on memory location.
     * The location is defined by how the matrix is paritioned among processes.
     * @tparam Tensor the backend type of the 2nd order tensor
     * @return the row range of the input matrix v.
     */
    template<class Derived>
    inline Range row_range(const Tensor<Derived, 2> &v) {
        return v.derived().row_range();
    }

    /*!
     * @fn Range col_range(const Tensor<Derived, 2> &v)
     * @brief Column ranges allow to coordinate the access of columns of the tensor based on memory location.
     * The location is defined by how the matrix is paritioned among processes.
     * @tparam Tensor the backend type of the 2nd order tensor
     * @return the column range of the input matrix v.
     */

    template<class Derived>
    inline Range col_range(const Tensor<Derived, 2> &v) {
        return v.derived().col_range();
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
    template<class Derived, int Order>
    inline bool read(const std::string &path, Tensor<Derived, Order> &t) {
        return t.derived().read(path);
    }

    /**
     * @ingroup    io
     * @brief      Writes and saves tensor into file.
     *
     * @param[in]  path    The path.
     * @param[in]  t       Tensor to be saved.
     */
    template<class Derived, int Order>
    inline bool write(const std::string &path, const Tensor<Derived, Order> &t) {
        return t.derived().write(path);
    }


    /**
     * @ingroup    interoperability
     * @brief      Converts utopia-wrapper specified type of tensor into original type.
     *
     * @param      t        Utopia tensor(wrapper).
     * @param      rawType  The external library tensor (backend supported).
     */
    template<class Derived, int Order, class RawType>
    inline void convert(Tensor<Derived, Order> &/*t*/, RawType &/*rawType*/) {
       assert(false && "OVERLOAD ME IF NEEDED");
    }

    /**
     * @ingroup    interoperability
     * @brief      Converts bewtween vectors with compatible types
     *
     * @param      from  source of copy
     * @param      to    destination of copy
     */
    template<typename T1, typename T2>
    inline void convert(const std::vector<T1> &from, std::vector<T2> &to) {
       to.resize(from.size());
       std::copy(from.begin(), from.end(), to.begin());
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
    template<class Derived, int Order>
    void disp(const Tensor<Derived, Order> &w)
    {
        w.derived().describe();
    }

    template<class T>
    void disp(const std::vector<T> &w)
    {
        for(const auto &wi : w) {
            disp(wi);
            disp("\n");
        }
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
    template<class Derived, int Order>
    bool empty(const Tensor<Derived, Order> &w)
    {
        return w.derived().empty();
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
    
    template<class Derived, int Order>
    bool has_nan_or_inf(const Tensor<Derived, Order> &w)
    {
        return w.derived().has_nan_or_inf();
    }

    template<class Derived>
    inline constexpr int order(const Expression<Derived> &)
    {
        return Derived::Order;
    }
}

#endif //utopia_utopia_MATRIXEXPR_HPP
