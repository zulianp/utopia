//
// Created by Patrick Zulian on 28/10/15.
//

#ifndef UTOPIA_UTOPIA_DIAG_HPP
#define UTOPIA_UTOPIA_DIAG_HPP

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include <string>

namespace utopia {
    template<class Expr>
    class Diag : public Expression<Diag<Expr> > {
    public:
        Diag(const Expr &expr) : expr_(expr)
        {}

        typedef typename Expr::Scalar Scalar;

        enum {
            Order = (Expr::Order == 1) ? 2 : 1
        };

        inline const Expr &expr() const { return expr_; }

        std::string getClass() const {
            return  "Diag<" + expr_.getClass() + ">";
        }

    private:
       UTOPIA_STORE_CONST(Expr) expr_;
        // const Expr &expr_;
    };

    template<class Expr>
    class Traits< Diag<Expr> > : public Traits<Expr> {};


    /**
     * @ingroup     reductions
     * @brief       \f$  diag( \cdot ) \f$ 
     */
    template<class Derived>
    Diag<Derived> diag(const Expression<Derived> &expr) {
        return Diag<Derived>(expr.derived());
    }

    template<class Expr>
    Size size(const Diag<Expr> &diag)
    {
        if(Expr::Order == 2) {
            Size s = size(diag.expr());
            return { std::min(s.get(0), s.get(1)) };
        } else {
            assert(Expr::Order == 1);
            Size s = size(diag.expr());
            return { std::min(s.get(0), s.get(0)) };
        }
    }
}

#endif //UTOPIA_UTOPIA_DIAG_HPP
