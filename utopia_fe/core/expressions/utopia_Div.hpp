#ifndef UTOPIA_FE_DIV_HPP
#define UTOPIA_FE_DIV_HPP

#include "utopia_DifferentialOperator.hpp"
#include "utopia_FunctionalTraits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Expr>
    class Divergence : public DifferentialOperator<Divergence<Expr> > {
    public:
        enum { Order = Expr::Order - 1 };

        typedef typename Expr::Scalar Scalar;

        std::string get_class() const override { return "Divergence<" + expr_.get_class() + ">"; }

        inline const Expr &expr() const { return expr_; }

        Divergence(const Expr &expr) : expr_(expr) {}

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template <class Derived>
    inline Divergence<Derived> div(const Expression<Derived> &expr) {
        return Divergence<Derived>(expr);
    }

    template <class Expr>
    class Traits<Divergence<Expr> > : public Traits<Expr> {
    public:
        enum { FILL_TYPE = utopia::FillType::DENSE };
    };

    template <class Expr, class AssemblyContext>
    class FunctionalTraits<Divergence<Expr>, AssemblyContext> {
    public:
        inline static int type(const Divergence<Expr> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);
        }
        inline static int order(const Divergence<Expr> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx) - 1;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_DIV_HPP
