//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_REDUCE_HPP
#define utopia_utopia_REDUCE_HPP

#include <string>
#include "utopia_Expression.hpp"
#include "utopia_Operators.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Size.hpp"

namespace utopia {

    template<class Expr_, class Operation>
    class Reduce : public Expression< Reduce<Expr_, Operation> > {
    public:
        typedef Expr_ Expr;
        typedef typename Expr::Scalar Scalar;

        enum {
            Order = 0
        };

        inline Reduce(const Expr &expr, const Operation &operation = Operation())
                : _expr(expr), _operation(operation) { }

        inline std::string getClass() const {
            return "Reduce<" + _expr.getClass() + ", " + _operation.getClass() + ">";
        }

        inline const Expr &expr() const
        {
            return _expr;
        }

        inline const Operation &operation() const
        {
            return _operation;
        }


       operator typename Traits<Reduce>::Scalar() const
       {
           Evaluator<typename Traits<Reduce>::Vector, Traits<Reduce>::Backend> eval;
           return eval.eval(*this);
       }


    private:
        UTOPIA_STORE_CONST(Expr) _expr;
        Operation _operation;
    };

    template<class Expr, class Operation>
    class Traits< Reduce<Expr, Operation> > : public Traits<Expr> {};


    template<class Derived, class Operation>
    inline Reduce<Derived, Operation> reduce(const Expression<Derived> &expr, const Operation &operation) {
        return Reduce<Derived, Operation>(expr.derived(), operation);
    }


     /**     @defgroup   reductions Reductions
     *       @ingroup    algebra
     */

     /**
     * @ingroup    reductions
     * @brief      \f$   \sum_{i = 0}^{n = dim(x)} x_i   \f$
     */
    template<class Derived>
    inline Reduce<Derived, Plus> sum(const Expression<Derived> &expr) {
        return Reduce<Derived, Plus>(expr.derived());
    }

    template<class Derived>
    inline Reduce<Derived, Min> min(const Expression<Derived> &expr) {
        return Reduce<Derived, Min>(expr.derived());
    }

    /**
     * @ingroup     reductions
     * @brief       \f$  <x, y> \f$
     */
    template<class Left, class Right>
    inline Reduce<Binary<Left, Right, EMultiplies>, Plus> dot(const Expression<Left> &left, const Expression<Right> &right) {
        return sum(Binary<Left, Right, EMultiplies>(left.derived(), right.derived()));
    }

    template<class Left, class Right>
    inline Reduce<Binary<Left, Right, EMultiplies>, Plus> inner(const Expression<Left> &left, const Expression<Right> &right) {
        return sum(Binary<Left, Right, EMultiplies>(left.derived(), right.derived()));
    }


//    template<class Derived>
//    inline Reduce<Derived, AbsPlus> norm1(const Expression<Derived> &expr) {
//
//        return Reduce<Derived, AbsPlus>(expr.derived());
//    }

    template<class Expr, class Operation>
    inline Size size(const Reduce<Expr, Operation> &expr)
    {
        Size s(1);
        s.set(0, 1);
        return s;
    }

    //Other functions: norm2(.) norm1(.), normInf(.) ...

}

#endif //utopia_utopia_REDUCE_HPP
