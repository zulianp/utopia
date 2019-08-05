#ifndef UTOPIA_DIFFERENTIABLE_HPP
#define UTOPIA_DIFFERENTIABLE_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<class Expr>
    class Differentiable : public Expression< Differentiable<Expr> > {
    public:
        typedef typename Expr::Scalar Scalar;

        static const int Order = Expr::Order;

        Differentiable(const Expr &expr)
        : expr_(expr)
        {}

        inline const Expr &expr() const
        {
            return expr_;
        }

        inline std::string getClass() const {
            return "independent_(" + expr_.getClass() + ")";
        }

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template<class Expr>
    inline Size size(const Differentiable<Expr> &expr)
    {
        return size(expr.expr());
    }

    template<class InnerExpr>
    class TreeProperties< Differentiable<InnerExpr> > {
    public:
        enum { greatest_tensor_order = TreeProperties<InnerExpr>::greatest_tensor_order };
        enum { smallest_tensor_order = TreeProperties<InnerExpr>::smallest_tensor_order };
        enum { has_mat_mat_mul 		 = TreeProperties<InnerExpr>::has_mat_mat_mul 		};
        enum { has_differentiable_sub_tree = 1 };
    };

    template<class Expr>
    class Traits< Differentiable<Expr> > : public Traits<Expr> {};

    template<class Expr>
    constexpr int is_differentiable(const Expr &expr)
    {
        return TreeProperties<Expr>::has_differentiable_sub_tree;
    }

    template<class Expr>
    constexpr int is_differentiable()
    {
        return TreeProperties<Expr>::has_differentiable_sub_tree;
    }

    template<class Derived>
    inline Differentiable<Derived> independent_variable(const Expression<Derived> &expr)
    {
        return expr.derived();
    }
}

#endif //UTOPIA_DIFFERENTIABLE_HPP

