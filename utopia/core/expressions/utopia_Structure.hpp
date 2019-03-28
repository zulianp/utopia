#ifndef UTOPIA_STRUCTURE_HPP
#define UTOPIA_STRUCTURE_HPP

#include "utopia_Expression.hpp"
#include "utopia_Wrapper.hpp"

namespace utopia {

    template<class Expr>
    class Structure final : public Expression<Structure<Expr>> {
    public:
        std::string getClass() const override {
            return "Structure<" + expr_.getClass() + ">";
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

    template<class Tensor>
    Structure<Wrapper<Tensor, 2>> structure(const Wrapper<Tensor, 2> &mat)
    {
        return Structure<Wrapper<Tensor, 2>>(mat);
    }

    template<class Tensor>
    Structure<Wrapper<Tensor, 1>> structure(const Wrapper<Tensor, 1> &mat)
    {
        return Structure<Wrapper<Tensor, 1>>(mat);
    }
}

#endif //UTOPIA_STRUCTURE_HPP
