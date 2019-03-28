/*
* @Author: alenakopanicakova
* @Date:   2016-08-016
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2016-08-16
*/

#ifndef UTOPIA_MATRIX_RESTRICTION_HPP
#define UTOPIA_MATRIX_RESTRICTION_HPP

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include <string>

namespace utopia {

    /**
     * @defgroup   tensor_products Tensor Products
     * @ingroup     algebra
     */

     /**
     * @ingroup     tensor_products
     * @brief       Creates product \f$ P^T * A * P \f$ - the same exact procedure and performance can be achieved by calling \f$ P^T * A * P \f$ directly. \n
     */
    template<class Left, class Right>
    Multiply< Multiply<Transposed<Right>, Left>, Right> ptap(const Expression<Left> &A, const Expression<Right> &P)
    {
        return transpose(P) * A * P;
    }
}

#endif //UTOPIA_MATRIX_RESTRICTION_HPP
