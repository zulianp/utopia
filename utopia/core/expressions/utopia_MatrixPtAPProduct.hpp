#ifndef UTOPIA_MATRIX_RESTRICTION_HPP
#define UTOPIA_MATRIX_RESTRICTION_HPP

#include <string>
#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"

namespace utopia {

    /**
     * @defgroup   tensor_products Tensor Products
     * @ingroup     algebra
     */

    /**
     * @ingroup     tensor_products
     * @brief       Creates product \f$ P^T * A * P \f$ - the same exact procedure and performance can be achieved by
     * calling \f$ P^T * A * P \f$ directly. \n
     */
    template <class Left, class Right>
    Multiply<Multiply<Transposed<Right>, Left>, Right> ptap(const Expression<Left> &A, const Expression<Right> &P) {
        return transpose(P) * A * P;
    }
}  // namespace utopia

#endif  // UTOPIA_MATRIX_RESTRICTION_HPP
