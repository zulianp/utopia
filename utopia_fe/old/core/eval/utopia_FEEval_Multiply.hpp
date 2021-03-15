#ifndef UTOPIA_FE_EVAL_MULTIPLY_HPP
#define UTOPIA_FE_EVAL_MULTIPLY_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
    template <class Left, class Function, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Left, Gradient<Function>>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<Left, Gradient<Function>> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
                                                     expr.right(),
                                                     ctx)) {
            return FEBackend<Backend>::multiply(
                FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx), expr.right(), ctx);
        }
    };

    template <class Left, class Function, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Left, TestFunction<Function>>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<Left, TestFunction<Function>> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
                                                     expr.right(),
                                                     ctx)) {
            return FEBackend<Backend>::multiply(
                FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx), expr.right(), ctx);
        }
    };

    template <class Left, class Function, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Left, TrialFunction<Function>>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<Left, TrialFunction<Function>> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
                                                     expr.right(),
                                                     ctx)) {
            return FEBackend<Backend>::multiply(
                FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx), expr.right(), ctx);
        }
    };

    template <class Left, class Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Left, Right>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Multiply<Left, Right> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
                                                     FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(),
                                                                                                       ctx),
                                                     ctx)) {
            return FEBackend<Backend>::multiply(FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
                                                FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx),
                                                ctx);
        }
    };

    // template<class C, class F, class Right, class Traits, int Backend, int IsQuadData>
    // class FEEval<Multiply<Interpolate<C, F>, Right>, Traits, Backend, IsQuadData> {
    // public:
    // 	typedef utopia::Interpolate<C, F> Left;
    // 	typedef utopia::Multiply<Left, Right> Expr;

    //     inline static auto apply(
    //     	const Expr &expr,
    //     	AssemblyContext<Backend> &ctx) -> decltype(
    //     	FEBackend<Backend>::multiply(
    //     		FEBackend<Backend>::fun(expr.left(), ctx),
    //     		FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx), ctx))
    //     {
    //     	return FEBackend<Backend>::multiply(
    //     		FEBackend<Backend>::fun(expr.left(), ctx),
    //     		FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx), ctx);
    //     }
    // };

    template <class C, class F, class Space, class Traits, int Backend, int IsQuadData>
    class FEEval<Multiply<Interpolate<C, F>, TrialFunction<Space>>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Interpolate<C, F> Left;
        typedef utopia::TrialFunction<Space> Right;
        typedef utopia::Multiply<Left, Right> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::multiply(FEBackend<Backend>::fun(expr.left(), ctx),
                                                     FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(),
                                                                                                       ctx),
                                                     ctx)) {
            return FEBackend<Backend>::multiply(FEBackend<Backend>::fun(expr.left(), ctx),
                                                FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx),
                                                ctx);
        }
    };

    template <class Left, class Right, class AssemblyContext>
    class FunctionalTraits<Multiply<Left, Right>, AssemblyContext> {
    public:
        inline static int type(const Multiply<Left, Right> &expr, const AssemblyContext &ctx) {
            return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx),
                            FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
        }

        inline static int order(const Multiply<Left, Right> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
                   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_MULTIPLY_HPP
