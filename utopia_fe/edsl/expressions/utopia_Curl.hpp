#ifndef UTOPIA_FE_CURL_HPP
#define UTOPIA_FE_CURL_HPP

#include "utopia_DifferentialOperator.hpp"
#include "utopia_FunctionalTraits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Expr>
    class Curl : public DifferentialOperator<Curl<Expr> > {
    public:
        enum { Order = Expr::Order - 1 };

        typedef typename Expr::Scalar Scalar;

        std::string get_class() const override { return "Curl<" + expr_.get_class() + ">"; }

        inline const Expr &expr() const { return expr_; }

        Curl(const Expr &expr) : expr_(expr) {}

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template <class Derived>
    inline Curl<Derived> curl(const Expression<Derived> &expr) {
        return Curl<Derived>(expr);
    }

    template <class Expr>
    class Traits<Curl<Expr> > : public Traits<Expr> {
    public:
        enum { FILL_TYPE = utopia::FillType::DENSE };
    };

    template <class Expr, class AssemblyContext>
    class FunctionalTraits<Curl<Expr>, AssemblyContext> {
    public:
        inline static int type(const Curl<Expr> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);
        }
        inline static int order(const Curl<Expr> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx) - 1;
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_CURL_HPP
