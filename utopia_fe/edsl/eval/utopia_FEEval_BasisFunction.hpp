#ifndef UTOPIA_FE_EVAL_BASIS_FUNCTION_HPP
#define UTOPIA_FE_EVAL_BASIS_FUNCTION_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEEval_Empty.hpp"
#include "utopia_fe_lang.hpp"

namespace utopia {

    template <class Derived, class FunctionSpaceT, class Traits, int Backend, int IsQuadData>
    class FEEval<BasisFunction<Derived, FunctionSpaceT>, Traits, Backend, IsQuadData> {
    public:
        typedef BasisFunction<Derived, FunctionSpaceT> Expr;
        typedef EXPR_TYPE(Traits, Expr) Result;

        inline static Result apply(const Expr &expr, AssemblyContext<Backend> &ctx) {
            assert(false);
            return Result();
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_BASIS_FUNCTION_HPP
