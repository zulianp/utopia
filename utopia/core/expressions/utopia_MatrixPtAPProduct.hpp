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

namespace utopia 
{
    template<class Left, class Right>
    class MatrixPtAPProduct : public Expression< MatrixPtAPProduct<Left, Right> > 
    {
    public:
        MatrixPtAPProduct(const Left &left, const Right &right)
                : _left(left), _right(right)
        {}

        inline const Left &left() const 
        {
            return _left;
        }

        inline const Right &right() const 
        {
            return _right;
        }

        std::string getClass() const override
        {
            return "MatrixPtAPProduct<" + left().getClass() + ", " + right().getClass() + ">";
        }


    private:
        UTOPIA_STORE_CONST(Left)  _left;
        UTOPIA_STORE_CONST(Right) _right;
    };

    /**
     * @defgroup   tensor_products Tensor Products
     * @ingroup     algebra
     */

     /**
     * @ingroup     tensor_products
     * @brief       Creates product \f$ P^T * A * P \f$ - the same exact procedure and performance can be achieved by calling \f$ P^T * A * P \f$ directly. \n
     */
    template<class Left, class Right>
    inline MatrixPtAPProduct<Left, Right> mat_PtAP_product(const Expression<Left> &left, const Expression<Right> &shape_vector)
    {
        return MatrixPtAPProduct<Left, Right>(left.derived(), shape_vector.derived());
    }
}

#endif //UTOPIA_MATRIX_RESTRICTION_HPP
