#ifndef UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP

#include <iostream>
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Reduce.hpp"
#include "utopia_Assign.hpp"
#include "utopia_Structured.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Ranged.hpp"
#include "utopia_Multiply.hpp"
#include "utopia_Transposed.hpp"
#include "utopia_Boolean.hpp"
#include "utopia_OuterProduct.hpp"
#include "utopia_Norm.hpp"
#include "utopia_FillTypeQuery.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Macros.hpp"

namespace utopia {

    template<class Expr, class Traits = utopia::Traits<Expr>, int Backend = Traits::Backend>
    class Eval {};

    template<class Derived, int Order, class Traits, int Backend>
    class Eval<Tensor<Derived, Order>, Traits, Backend> {
    public:
        inline static const Derived &apply(const Tensor<Derived, Order> &expr) {
            return expr.derived();
        }

        inline static Derived &apply(Tensor<Derived, Order> &expr) {
            return expr.derived();
        }

        inline static Derived && apply(Tensor<Derived, Order> &&expr) {
            return std::move(expr.derived());
        }

        inline static void apply(Tensor<Derived, Order> &&expr, Derived &result) {
            result = std::move(expr.derived());
        }

        inline static void apply(const Tensor<Derived, Order> &expr, Derived &result) {
            result = expr.derived();
        }
    };

    template<class Derived, int Order, class Traits, int Backend>
    class Eval<Tensor<Derived &, Order>, Traits, Backend> {
    public:
        inline static const Derived &apply(const Tensor<Derived &, Order> &expr) {
            return expr.derived();
        }

        inline static Derived &apply(Tensor<Derived &, Order> &expr) {
            return expr.derived();
        }
    };

    template<class Derived, int Order, class Traits, int Backend>
    class Eval<Tensor<const Derived &, Order>, Traits, Backend> {
    public:
        inline static const Derived &apply(const Tensor<const Derived &, Order> &expr) {
            return expr.derived();
        }
    };

    //for c++14 only
    // template<class Expr>
    // auto eval(const Expression<Expr> &expr)
    // {
    //     return Eval<Expr>::apply(expr.derived());
    // }
}

#endif //UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP
