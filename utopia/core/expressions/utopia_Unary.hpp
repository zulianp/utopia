//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef SIMMOD_UNARY_HPP
#define SIMMOD_UNARY_HPP

#include "utopia_Expression.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template<class _Expr, class _Operation>
    class Unary : public Expression< Unary<_Expr, _Operation> > {
    public:
        using Expr = _Expr;
        using Operation = _Operation;
        using Scalar = typename utopia::Traits<Expr>::Scalar;

        static const int Order = Expr::Order;

        Unary(const Expr &expr, const Operation operation = Operation()) : _expr(expr), _operation(operation) { }
        inline const Expr &expr() const { return _expr; }

        virtual ~Unary() = default;

        std::string get_class() const override { return _operation.get_class() + "<" + _expr.get_class() + ">"; }

        inline const Operation &operation() const
        {
            return _operation;
        }

    private:
        UTOPIA_STORE_CONST(Expr) _expr;
        Operation _operation;
    };

    template<class Expr, class Operation>
    class Traits< Unary<Expr, Operation> > : public Traits<Expr> {};


    template<class Expr, class Operation>
    auto size(const Unary<Expr, Operation> &expr) -> decltype( size(expr.expr()) ) {
        return size(expr.expr());
    }
}

#endif //SIMMOD_UNARY_HPP
