#ifndef UTOPIA_INVERSE_HPP
#define UTOPIA_INVERSE_HPP

#include "utopia_Expression.hpp"

#include <type_traits>

namespace utopia {

    template <class Expr_>
    class Inverse : public Expression<Inverse<Expr_> > {
    public:
        using Expr = Expr_;
        using Scalar = typename std::remove_const<typename Expr::Scalar>::type;

        static const int Order = Expr::Order;

        // static_assert(Order == 2, "Only works for matrices");

        Inverse(const Expr &expr) : expr_(expr) {}

        inline const Expr &expr() const { return expr_; }

        std::string get_class() const override { return "Inverse<" + expr_.get_class() + ">"; }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template <class Expr>
    class Traits<Inverse<Expr> > : public Traits<Expr> {};

    template <class Derived>
    Inverse<Derived> inv(const Expression<Derived> &expr) {
        return Inverse<Derived>(expr.derived());
    }

    template <class Expr>
    Size size(const Inverse<Expr> &expr) {
        auto s = size(expr.expr());
        return s;
    }
}  // namespace utopia

#endif  // UTOPIA_INVERSE_HPP
