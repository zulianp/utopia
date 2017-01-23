//
// Created by Patrick Zulian on 28/05/15.
//

#ifndef UTOPIA_UTOPIA_TRANSPOSE_HPP
#define UTOPIA_UTOPIA_TRANSPOSE_HPP

#include "utopia_Expression.hpp"

namespace utopia {
    template<class _Expr>
    class Transposed : public Expression< Transposed<_Expr> > {
    public:
        typedef _Expr Expr;
        typedef typename Expr::Scalar Scalar;

        enum {
            Order = Expr::Order
        };

        Transposed(const Expr &expr) : _expr(expr) { }
        inline const Expr &expr() const { return _expr; }

        virtual ~Transposed() { }

        std::string getClass() const {
            return "Transposed<" + _expr.getClass() + ">";
        }

    private:
        UTOPIA_STORE_CONST(Expr) _expr;
    };

    template<class Expr>
    class Traits< Transposed<Expr> > : public Traits<Expr> {};

    /**     @defgroup   permutation Permutations
     *      @ingroup    algebra
     */

    /**
     * @ingroup     permutation
     * @brief       \f$ A^T \f$
     */
    template<class Derived>
    Transposed<Derived> transpose(const Expression<Derived> &expr) {
        return Transposed<Derived>(expr.derived());
    }

    template<class Derived>
    bool is_transposed(const Expression<Derived> &)
    {
        return false;
    }

    template<class Expr>
    bool is_transposed(const Transposed<Expr> &)
    {
        return true;
    }

    template<class Expr>
    Size size(const Transposed<Expr> &expr)
    {
        auto s = size(expr.expr());
        if(Expr::Order == 1) {
            return s;
        }

        Size ret(2);
        ret.set(0, s.get(1));
        ret.set(1, s.get(0));
        return ret;
    }


}

#endif //UTOPIA_UTOPIA_TRANSPOSE_HPP
