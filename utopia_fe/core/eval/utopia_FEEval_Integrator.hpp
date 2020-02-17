#ifndef UTOPIA_FE_EVAL_INTEGRATOR_HPP
#define UTOPIA_FE_EVAL_INTEGRATOR_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_Integrator.hpp"

namespace utopia {

    template<class FunctionSpace, class Traits, int Backend, int IsQuadData>
    class FEEval< LinearIntegrator<FunctionSpace>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::LinearIntegrator<FunctionSpace> Expr;
        using Type = typename Expr::Type;

        inline static Type apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) 
        {
            Type result;
            expr.assemble(ctx, result);
            return result;
        }
    };

    template<class FunctionSpace, class Traits, int Backend, int IsQuadData>
    class FEEval< BilinearIntegrator<FunctionSpace>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::BilinearIntegrator<FunctionSpace> Expr;
        using Type = typename Expr::Type;

        inline static Type apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) 
        {
            Type result;
            expr.assemble(ctx, result);
            return result;
        }
    };

}

#endif //UTOPIA_FE_EVAL_INTEGRATOR_HPP
