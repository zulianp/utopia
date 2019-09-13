//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP

#include <iostream>
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Reduce.hpp"
#include "utopia_Backend.hpp"
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
// #include "utopia_Tracer.hpp"
// #include "utopia_Tensor.hpp"


namespace utopia {

#define UTOPIA_BACKEND(Traits_) (utopia::Backend<typename Traits_::Scalar, Traits_::Backend>::Instance())

    template<class Expr, class Traits = utopia::Traits<Expr>, int Backend = Traits::Backend>
    class Eval {};

    // template<class Tensor, int Order, class Traits, int Backend>
    // class Eval<Wrapper<Tensor, Order>, Traits, Backend> {
    // public:
    //     inline static const Tensor &apply(const Wrapper<Tensor, Order> &expr) {
    //         return expr.implementation();
    //     }

    //     inline static Tensor &apply(Wrapper<Tensor, Order> &expr) {
    //         return expr.implementation();
    //     }
    // };

    // template<class Tensor, int Order, class Traits, int Backend>
    // class Eval<Wrapper<Tensor &, Order>, Traits, Backend> {
    // public:
    //     inline static const Tensor &apply(const Wrapper<Tensor &, Order> &expr) {
    //         return expr.implementation();
    //     }

    //     inline static Tensor &apply(Wrapper<Tensor &, Order> &expr) {
    //         return expr.implementation();
    //     }
    // };

    // template<class Tensor, int Order, class Traits, int Backend>
    // class Eval<Wrapper<const Tensor &, Order>, Traits, Backend> {
    // public:
    //     inline static const Tensor &apply(const Wrapper<const Tensor &, Order> &expr) {
    //         return expr.implementation();
    //     }
    // };


    template<class Derived, int Order, class Traits, int Backend>
    class Eval<Tensor<Derived, Order>, Traits, Backend> {
    public:
        inline static const Derived &apply(const Tensor<Derived, Order> &expr) {
            return expr.derived();
        }

        inline static Derived &apply(Tensor<Derived, Order> &expr) {
            return expr.derived();
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


}

#endif //UTOPIA_UTOPIA_EVAL_EMPTY_HPP_HPP
