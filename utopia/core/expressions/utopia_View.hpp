#ifndef UTOPIA_UTOPIA_VIEW_HPP
#define UTOPIA_UTOPIA_VIEW_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Mutable.hpp"
#include "utopia_Range.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Expr>
    class View : public Expression<View<Expr> > {
    public:
        static const int Order = Expr::Order;
        using Scalar = typename utopia::Traits<Expr>::Scalar;

        template <class Derived>
        View &operator=(const Expression<Derived> &expr) {
            using A = utopia::Assign<View, Derived>;

            EvalAssignToView<Expr, Derived, Traits<Expr>, Traits<Expr>::Backend>::apply(A(*this, expr.derived()));

            return *this;
        }

        View(Expr &expr, const Range &row_range, const Range &col_range)
            : _expr(expr), row_range_(row_range), col_range_(col_range) {}

        inline Expr &expr() const { return _expr; }

        friend const Range &row_range(const View &v) { return v.row_range_; }

        friend const Range &col_range(const View &v) { return v.col_range_; }

    private:
        Expr &_expr;
        const Range row_range_;
        const Range col_range_;
    };

    template <class Expr>
    class Traits<View<Expr> > : public Traits<Expr> {};
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_VIEW_HPP