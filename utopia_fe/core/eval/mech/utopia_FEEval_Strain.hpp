#ifndef UTOPIA_FE_EVAL_STRAIN_HPP
#define UTOPIA_FE_EVAL_STRAIN_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Interpolate.hpp"

namespace utopia {

    template<class C, class LeftFun, class RightFun>
    using StrainExpr = 
        utopia::Binary<
            Transposed<
                GradInterpolate<C, LeftFun>
                >,
            GradInterpolate<C, RightFun>,
            Plus
        >;
    

    template<class C, class LeftFun, class RightFun, class Traits, int Backend, int IsQuadData>
    class FEEval< StrainExpr<C, LeftFun, RightFun>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::StrainExpr<C, LeftFun, RightFun> Expr;
        typedef utopia::BinaryDelegate<Expr, Traits, IsQuadData> DelegateT;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEBackend<Backend>::grad_t_plus_grad(expr.left().expr(), expr.right(), ctx) )
        {
            return FEBackend<Backend>::grad_t_plus_grad(
                1.0,
                expr.left().expr(),
                expr.right(),
                ctx
            );
        }
    };
}

#endif //UTOPIA_FE_EVAL_STRAIN_HPP
