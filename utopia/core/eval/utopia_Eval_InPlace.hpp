//
// Created by Patrick Zulian on 30/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_INPLACE_HPP
#define UTOPIA_UTOPIA_EVAL_INPLACE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {
    template<class Left, class Right, class Operation, class Traits, int Backend>
    class Eval<InPlace<Left, Right, Operation>, Traits, Backend> {
    public:

        inline static bool apply(const InPlace<Left, Right, Operation> &expr)
        {
            //FIXME connect to backend without tree transformation
            typedef utopia::Binary<Left, Right, Operation> TransformedExpr;
            typedef utopia::Assign<Left, TransformedExpr> AssignExpr;

            return Eval<AssignExpr, Traits>::apply(
                    AssignExpr(
                            expr.left(),
                            TransformedExpr(
                                    expr.left(),
                                    expr.right(),
                                    expr.operation()
                            )
                    )
            );
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_INPLACE_HPP
