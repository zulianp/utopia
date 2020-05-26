#ifndef UTOPIA_FE_EVAL_EMPTY_HPP
#define UTOPIA_FE_EVAL_EMPTY_HPP

#include "utopia_Eval.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {
    // typedef int EXPR_CTX;
    static const int QUAD_DATA_YES = 1;
    static const int QUAD_DATA_NO = 0;

    template <class Expr, class Traits, int Backend, int IsQuadData>
    class FEEval : public Eval<Expr, Traits, Backend> {
    public:
        // default fallback on eval if does not exists
        inline static auto apply(const Expr &expr, const AssemblyContext<Backend> &)
            -> decltype(Eval<Expr, Traits, Backend>::apply(expr)) {
            return Eval<Expr>::apply(expr);
        }
    };

    template <class Derived, int Order, class Traits, int Backend, int IsQuadData>
    class FEEval<Tensor<Derived, Order>, Traits, Backend, IsQuadData> {
    public:
        // default fallback on eval if does not exists
        inline static const Derived &apply(const Tensor<Derived, Order> &expr, const AssemblyContext<Backend> &) {
            return expr.derived();
        }

        inline static Derived &apply(Tensor<Derived, Order> &expr, const AssemblyContext<Backend> &) {
            return expr.derived();
        }

        inline static Derived apply(Tensor<Derived, Order> &&expr, const AssemblyContext<Backend> &) {
            return std::move(expr.derived());
        }
    };

    // template<class Expr, class Traits, int Backend, int IsQuadData>
    // class NonTerminalFEEval : public FEEval<Expr, Traits, Backend, IsQuadData> {};

}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_EMPTY_HPP
