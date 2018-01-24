//
// Created by Patrick Zulian on 06/07/15.
//

#ifndef UTOPIA_UTOPIA_OUTERPRODUCT_HPP
#define UTOPIA_UTOPIA_OUTERPRODUCT_HPP

#include "utopia_StoreAs.hpp"
#include "utopia_Expression.hpp"

namespace utopia {
    template<class Left, class Right>
    class OuterProduct : public Expression< OuterProduct<Left, Right> > {
    public:
        static const int Order =  Left::Order + Right::Order;

        enum {
            StoreAs = UTOPIA_BY_VALUE
        };


        OuterProduct(const Left &left, const Right &right)
                : _left(left), _right(right)
        {}

        inline const Left &left() const {
            return _left;
        }

        inline const Right &right() const {
            return _right;
        }

        std::string getClass() const override
        {
            return "OuterProduct<" + left().getClass() + ", " + right().getClass() + ">";
        }

    private:
        UTOPIA_STORE_CONST(Left)  _left;
        UTOPIA_STORE_CONST(Right) _right;
    };

    /**
     * @ingroup     tensor_products
     * @brief       \f$  v_1 v_2^T  \f$
     */
    template<class Left, class Right>
    inline OuterProduct<Left, Right> outer(const Expression<Left> &left, const Expression<Right> &right)
    {
        return OuterProduct<Left, Right>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    class Traits< OuterProduct<Left, Right> > : public Traits<typename ChooseType<Left, Right, Right>::Type > {
    public:
        enum {
            FILL_TYPE = FillType::DENSE
        };
    };


    template<class Left, class Right>
    inline Size size(const OuterProduct<Left, Right> &expr)
    {
        Size result(2);
        result.set(0, size(expr.left()).get(0));
        result.set(1, size(expr.right()).get(0));
        return result;
    }
}

#endif //UTOPIA_UTOPIA_OUTERPRODUCT_HPP
