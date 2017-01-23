//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_VIEW_HPP
#define UTOPIA_UTOPIA_VIEW_HPP

#include "utopia_Mutable.hpp"
#include "utopia_Range.hpp"

namespace utopia {
    template<class Expr>
    class View : public Expression<View<Expr> >,
                 public Mutable<typename Expr::Implementation, View<Expr> > {
    public:
        enum {
            Order = Expr::Order
        };

        template<class Derived>
        View &operator=(const Expression<Derived> &expr) {
            _expr.evaluator().eval(Assign<View, Derived>(*this, expr.derived()));
            return *this;
        }

        View(Expr &expr, const Range &rowRange, const Range &colRange)
                : _expr(expr), _rowRange(rowRange), _colRange(colRange) { }

        Expr &expr() const {
            return _expr;
        }



        const Range &rowRange() const {
            return _rowRange;
        }

        const Range &colRange() const {
            return _colRange;
        }

        friend const Range &row_range(const View &v)
        {
            return v._rowRange;
        }

        friend const Range &col_range(const View &v)
        {
            return v._colRange;
        }

    private:
        Expr &_expr;
        const Range _rowRange;
        const Range _colRange;

    };

    template<class Expr>
    class Traits< View<Expr> > : public Traits<Expr> {};
}

#endif //UTOPIA_UTOPIA_VIEW_HPP