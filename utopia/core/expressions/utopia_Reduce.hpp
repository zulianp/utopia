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

        static const int Order = 0;

        inline Reduce(const Expr &expr, const Operation &operation = Operation())
                : expr_(expr), operation_(operation) { }

        inline std::string get_class() const override {
            return "Reduce<" + expr_.get_class() + ", " + operation_.get_class() + ">";
        }

        inline const Expr &expr() const
        {
            return expr_;
        }

        inline const Operation &operation() const
        {
            return operation_;
        }

       operator typename Traits<Reduce>::Scalar() const
       {
           return Eval<Reduce, Traits<Reduce>, Traits<Reduce>::Backend>::apply(*this);
       }


    private:
        UTOPIA_STORE_CONST(Expr) expr_;
        Operation operation_;
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

    template<class Derived>
    inline Reduce<Derived, Max> max(const Expression<Derived> &expr) {
        return Reduce<Derived, Max>(expr.derived());
    }




    template<class L, class R>
    using Dot = utopia::Reduce<Binary<L, R, utopia::EMultiplies>, utopia::Plus>;

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



    template<class Derived, typename T = typename Traits<Derived>::Scalar>
    inline Reduce<Derived, PlusIsNonZero<T>> nnz(const Expression<Derived> &expr, const T tol = 0.) {
        return Reduce<Derived, PlusIsNonZero<T>>(expr.derived(), PlusIsNonZero<T>(tol));
    }


//    template<class Derived>
//    inline Reduce<Derived, AbsPlus> norm1(const Expression<Derived> &expr) {
//
//        return Reduce<Derived, AbsPlus>(expr.derived());
//    }

    template<class Expr, class Operation>
    inline Size size(const Reduce<Expr, Operation> &/*expr*/)
    {
        Size s(1);
        s.set(0, 1);
        return s;
    }

    //Other functions: norm2(.) norm1(.), normInf(.) ...

}

#endif //utopia_utopia_REDUCE_HPP
