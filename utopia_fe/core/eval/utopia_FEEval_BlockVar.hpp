#ifndef UTOPIA_FE_EVAL_BLOCK_VAR_HPP
#define UTOPIA_FE_EVAL_BLOCK_VAR_HPP

#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
    template <typename T, class Traits, int Backend, int IsQuadData>
    class FEEval<BlockVar<T>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::BlockVar<T> Expr;

        inline static T apply(const Expr &expr, const AssemblyContext<Backend> &ctx) {
            return FEBackend<Backend>::apply(expr, ctx);
        }
    };

    template <typename Left, typename Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Binary<Number<Left>, BlockVar<Right>, Multiplies>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Binary<Number<Left>, BlockVar<Right>, Multiplies> Expr;

        inline static Right apply(const Expr &expr, const AssemblyContext<Backend> &ctx) {
            return static_cast<Left>(expr.left()) * FEBackend<Backend>::apply(expr.right(), ctx);
        }
    };

    template <typename Left, typename Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Binary<Number<Left>, BlockVar<Right>, Minus>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Binary<Number<Left>, BlockVar<Right>, Minus> Expr;

        inline static Right apply(const Expr &expr, const AssemblyContext<Backend> &ctx) {
            return static_cast<Left>(expr.left()) - FEBackend<Backend>::apply(expr.right(), ctx);
        }
    };

    template <typename Left, typename Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Binary<Number<Left>, BlockVar<Right>, Plus>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Binary<Number<Left>, BlockVar<Right>, Plus> Expr;

        inline static Right apply(const Expr &expr, const AssemblyContext<Backend> &ctx) {
            return static_cast<Left>(expr.left()) + FEBackend<Backend>::apply(expr.right(), ctx);
        }
    };

    template <typename Left, typename Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Binary<Number<Left>, BlockVar<Right>, Divides>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Binary<Number<Left>, BlockVar<Right>, Divides> Expr;

        inline static Right apply(const Expr &expr, const AssemblyContext<Backend> &ctx) {
            return static_cast<Left>(expr.left()) / FEBackend<Backend>::apply(expr.right(), ctx);
        }
    };

    template <class T, class AssemblyContext>
    class FunctionalTraits<BlockVar<T>, AssemblyContext> {
    public:
        inline static int type(const BlockVar<T> &expr, const AssemblyContext &) { return CONSTANT_FUNCTION; }

        inline static int order(const BlockVar<T> &expr, const AssemblyContext &) { return 0; }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_BLOCK_VAR_HPP
