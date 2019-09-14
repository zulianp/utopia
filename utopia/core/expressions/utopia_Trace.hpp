#ifndef UTOPIA_TRACE_HPP
#define UTOPIA_TRACE_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template<class Expr>
    class Trace : public Expression< Trace<Expr> > {
    public:
        static_assert(Expr::Order >= 2, "must be a 2nd order tensor or greater");

        static const int Order = 0;

        inline const Expr &expr() const
        {
            return expr_;
        }

        operator typename Traits<Trace>::Scalar() const
        {
            Evaluator<typename Traits<Trace>::Vector, Traits<Trace>::Backend> eval;
            return eval.eval(*this);
        }

        inline std::string get_class() const override
        {
            return "Trace<" + expr().get_class() + ">";
        }

        Trace(const Expr &expr) : expr_(expr) {}

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };


    template<class Expr>
    class Traits< Trace<Expr> > : public Traits<Expr> {};

    template<class Expr>
    inline Size size(const Trace<Expr> &/*expr*/)
    {
        return {1};
    }

    template<class Expr>
    class Unfold< Trace<Expr> > {
    public:
        typedef utopia::Reduce< utopia::Diag<Expr>,
                                utopia::Plus >  Type;

        inline static Type apply(const Trace<Expr> &expr)
        {
            return sum(diag(expr.expr()));
        }
    };

    template<class Expr>
    inline typename Unfold< Trace<Expr> >::Type shallow_unfold(const Trace<Expr> &expr)
    {
        return Unfold< Trace<Expr> >::apply(expr);
    }

     /**
     * @ingroup     reductions
     * @brief       \f$ tr(A) = \sum_{i = 0}^{n = dim(A)} A_{ii} \f$. \n
     *
     */
    template<class Derived>
    Trace<Derived> trace(const Expression<Derived> &expr)
    {
        static_assert(Derived::Order == 2, "expr must be a tensor of order 2");
        return expr.derived();
    }
}

#endif  //UTOPIA_TRACE_HPP
