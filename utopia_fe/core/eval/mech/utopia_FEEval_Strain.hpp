#ifndef UTOPIA_FE_EVAL_STRAIN_HPP
#define UTOPIA_FE_EVAL_STRAIN_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Interpolate.hpp"
#include "utopia_block_var.hpp"

namespace utopia {

    template <class C, class LeftFun, class RightFun>
    using StrainExpr = utopia::Binary<Transposed<GradInterpolate<C, LeftFun>>, GradInterpolate<C, RightFun>, Plus>;

    template <class C, class LeftFun, class RightFun, class Traits, int Backend, int IsQuadData>
    class FEEval<StrainExpr<C, LeftFun, RightFun>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::StrainExpr<C, LeftFun, RightFun> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::grad_t_plus_grad(1.0, expr.left().expr(), expr.right(), ctx)) {
            return FEBackend<Backend>::grad_t_plus_grad(1.0, expr.left().expr(), expr.right(), ctx);
        }
    };

    template <typename T, class C, class LeftFun, class RightFun, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<BlockVar<T>, StrainExpr<C, LeftFun, RightFun>>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<BlockVar<T>, StrainExpr<C, LeftFun, RightFun>> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::grad_t_plus_grad(FEBackend<Backend>::apply(expr.left(), ctx),
                                                             expr.right().left().expr(),
                                                             expr.right().right(),
                                                             ctx)) {
            return FEBackend<Backend>::grad_t_plus_grad(
                FEBackend<Backend>::apply(expr.left(), ctx), expr.right().left().expr(), expr.right().right(), ctx);
        }
    };

    ////////////////////////////////////////////////////////

    template <class C, class LeftFun, class RightFun, class Traits, int Backend, int IsQuadData>
    class FEEval<
        Multiply<Multiply<BlockVar<double>, Trace<StrainExpr<C, LeftFun, RightFun>>>, SymbolicTensor<Identity, 2>>,
        Traits,
        Backend,
        IsQuadData> {
    public:
        typedef utopia::Multiply<Multiply<BlockVar<double>, Trace<StrainExpr<C, LeftFun, RightFun>>>,
                                 SymbolicTensor<Identity, 2>>
            Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::grad_t_plus_grad(1.,
                                                             expr.left().right().expr().left().expr(),
                                                             expr.left().right().expr().right(),
                                                             ctx)) {
            auto &&l = FEBackend<Backend>::apply(expr.left().left(), ctx);
            auto &&r = expr.left().right().expr();

            return FEBackend<Backend>::trace_times_identity(
                FEBackend<Backend>::grad_t_plus_grad(l, r.left().expr(), r.right(), ctx), ctx);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_STRAIN_HPP
