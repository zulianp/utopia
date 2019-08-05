//
// Created by Patrick Zulian on 28/10/15.
//

#ifndef UTOPIA_UTOPIA_TENSOR_REDUCE_HPP
#define UTOPIA_UTOPIA_TENSOR_REDUCE_HPP

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include <string>
#include "utopia_Operators.hpp"

namespace utopia {
    template<class Expr, class Operation>
    class TensorReduce : public Expression< TensorReduce<Expr, Operation> > {
    public:
        TensorReduce(const Expr &expr, const int dim, const Operation op = Operation()) : expr_(expr), dim_(dim), op_(op)
        {}

        typedef typename Expr::Scalar Scalar;

        static_assert(Expr::Order == 2, "only supported for matrices");

        static const int Order = Expr::Order - 1;

        inline const Expr &expr() const { return expr_; }

        std::string getClass() const {
            return  "TensorReduce<" + expr_.getClass() + ">";
        }

        inline int dim() const
        {
            return dim_;
        }

        inline const Operation &operation() const
        {
            return op_;
        }

    private:
       UTOPIA_STORE_CONST(Expr) expr_;
       int dim_;
       Operation op_;
    };

    template<class Expr, class Operation>
    class Traits< TensorReduce<Expr, Operation> > : public Traits<Expr> {};


    template<class Derived>
    TensorReduce<Derived, Plus> sum(const Expression<Derived> &expr, const int dim) {
        return TensorReduce<Derived, Plus>(expr.derived(), dim);
    }

    template<class Derived>
    TensorReduce<Derived, Min> min(const Expression<Derived> &expr, const int dim) {
        return TensorReduce<Derived, Min>(expr.derived(), dim);
    }

    template<class Derived>
    TensorReduce<Derived, Max> max(const Expression<Derived> &expr, const int dim) {
        return TensorReduce<Derived, Max>(expr.derived(), dim);
    }


    template<class Expr, class Operation>
    inline Size size(const TensorReduce<Expr, Operation> &expr)
    {
        static_assert(Expr::Order == 2, "only supported for matrices");
        Size s = size(expr.expr());
        return { s.get(!expr.dim()) };
    }
}

#endif //UTOPIA_UTOPIA_TENSOR_REDUCE_HPP
