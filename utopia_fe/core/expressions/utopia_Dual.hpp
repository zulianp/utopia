#ifndef UTOPIA_DUAL_HPP
#define UTOPIA_DUAL_HPP

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FunctionalTraits.hpp"

namespace utopia {

    template<class Expr>
    class Dual : public Expression< Dual<Expr> > {
    public:
        static const int Order = Expr::Order;
        typedef typename Expr::Scalar Scalar;

        std::string get_class() const override { return "Dual<" + expr_.get_class() + ">"; }

        inline const Expr &expr() const
        {
            return expr_;
        }

        Dual(const Expr &expr)
        : expr_(expr)
        {}

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template<class Derived>
    inline Dual<Derived> dual(const Expression<Derived> &expr) {
        return Dual<Derived>(expr.derived());
    }

    template<class Expr>
    class Traits< Dual<Expr> > : public Traits<Expr> {
    public:
        enum {
            FILL_TYPE = utopia::FillType::DENSE
        };
    };

    template<class Expr, class AssemblyContext>
    class FunctionalTraits<Dual<Expr>, AssemblyContext>  {
    public:
        inline static int type(const Dual<Expr> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);  }
        inline static int order(const Dual<Expr> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx); }
    };


}

#endif //UTOPIA_DUAL_HPP
