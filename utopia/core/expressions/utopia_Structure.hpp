#ifndef UTOPIA_STRUCTURE_HPP
#define UTOPIA_STRUCTURE_HPP

#include "utopia_Expression.hpp"
#include "utopia_Tensor.hpp"

namespace utopia {

    template<class Expr>
    class Structure final : public Expression<Structure<Expr>> {
    public:
        std::string get_class() const override {
            return "Structure<" + expr_.get_class() + ">";
        }

        Structure(const Expr &expr) : expr_(expr) { }

        const Expr &expr() const
        {
            return expr_;
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template<class Expr>
    class Traits< Structure<Expr> > : public Traits<Expr> {};

    template<class Expr>
    auto size(const Structure<Expr> &expr) -> decltype( size(expr.expr()) ) {
        return size(expr.expr());
    }

    template<class Derived>
    Structure<Tensor<Derived, 2>> structure(const Tensor<Derived, 2> &mat)
    {
        return Structure<Tensor<Derived, 2>>(mat);
    }

    template<class Derived>
    Structure<Tensor<Derived, 1>> structure(const Tensor<Derived, 1> &mat)
    {
        return Structure<Tensor<Derived, 1>>(mat);
    }
}

#endif //UTOPIA_STRUCTURE_HPP
