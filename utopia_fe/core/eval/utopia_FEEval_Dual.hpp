#ifndef UTOPIA_FE_EVAL_DUAL_HPP
#define UTOPIA_FE_EVAL_DUAL_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Dual.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FEEval_Empty.hpp"

namespace utopia {

    template <class Tensor, class Traits, int Backend, int IsQuadData>
    class FEEval<Dual<Tensor>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Dual<Tensor> Expr;

        template <template <class> class Function, class Space>
        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::dual(expr.expr(), ctx)) {
            return FEBackend<Backend>::dual(expr.expr(), ctx);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_DUAL_HPP
