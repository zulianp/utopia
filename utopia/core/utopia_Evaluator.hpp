//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_EVALUATOR_HPP
#define utopia_utopia_EVALUATOR_HPP

#include <vector>

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


#include "utopia_Eval.hpp"

namespace utopia {
    template<class Result, int BAKEND_FLAG = Traits<Result>::Backend>
    class Evaluator {
    public:
        typedef utopia::Traits<Result> Traits;
        typedef typename Traits::Scalar  Scalar;

        typedef utopia::Backend<Scalar, Traits::Backend> BackendT;

        template<class Expr, int Order>
        static void eval(const Wrapper<Expr, Order> &expr, Size &size)
        {
            BackendT::Instance().size(eval(expr), size);
        }

        template<class Derived>
        static auto eval(const Expression<Derived> &expr) -> decltype(Eval<Derived, Traits>::apply(expr.derived()))
        {
            return Eval<Derived, Traits>::apply(expr.derived());
        }

        template<class Left, class Right>
        static void eval(const Construct<Left, Right> &expr)
        {
            Eval<Construct<Left, Right>, Traits>::apply(expr);
        }

        template<class Left, class Right>
        static void eval(const Assign<Left, Right> &expr)
        {
            Eval<Assign<Left, Right>, Traits>::apply(expr);
        }

        template<class Left, class Right, class Operation>
        static void eval(const InPlace<Left, Right, Operation> &expr)
        {
            Eval<InPlace<Left, Right, Operation>, Traits>::apply(expr);
        }      
    };
}

#endif //utopia_utopia_EVALUATOR_HPP
