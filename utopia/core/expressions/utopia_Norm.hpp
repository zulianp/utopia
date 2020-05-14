//
// Created by Patrick Zulian on 21/05/15.
//

#ifndef UTOPIA_UTOPIA_NORM_HPP
#define UTOPIA_UTOPIA_NORM_HPP
#include "utopia_Expression.hpp"
#include "utopia_Unfold.hpp"

namespace utopia {

    static const int INFINITY_NORM_TAG = -1;

    template<class Expr, int Type>
    class Norm : public Expression<Norm<Expr, Type> > {
    public:
        typedef typename Expr::Scalar Scalar;
        static const int Order = 0;

        Norm(const Expr &expr)
                : _expr(expr)
        {}

        inline const Expr &expr() const
        {
            return _expr;
        }

        std::string get_class() const override { return "Norm<" + _expr.get_class() + ">"; }

        operator typename Traits<Norm>::Scalar() const
        {
            // Evaluator<typename Traits<Norm>::Vector, Traits<Norm>::Backend> e;
            // return e.eval(*this);
            return Eval<Norm, Traits<Norm>, Traits<Norm>::Backend>::apply(*this);
        }

    private:
        UTOPIA_STORE_CONST(Expr) _expr;
    };

    template<class Expr>
    class Unfold< Norm<Expr, 2> > {
    public:
        typedef utopia::Unary<
                        utopia::Reduce<
                                utopia::Unary<Expr, utopia::Pow2>,
                                utopia::Plus >,
                        utopia::Sqrt
                        > Type;

        inline static Type apply(const Norm<Expr, 2> &expr)
        {
            return Type( Reduce <
                                Unary<Expr, Pow2>,
                                Plus
                                >( Unary<Expr, Pow2>(expr.expr()) )
                        );
        }
    };

    template<class Expr>
    class Unfold< Norm<Expr, 1> > {
    public:
        typedef utopia::Reduce<Expr, utopia::AbsPlus> Type;
    };

    // template<class Expr, int Order>
    // class Fold< typename Unfold< Norm<Expr, Order> >::Type > {
    // public:
    //     typedef utopia::Norm<Expr, Order> Type;

    //     inline static Type apply(const Norm<Expr, 2> &expr)
    //     {
    //         return expr.expr().expr().expr();
    //     }

    // };

    template<class Expr>
    inline typename Unfold< Norm<Expr, 2> >::Type shallow_unfold(const Norm<Expr, 2> &expr)
    {
        return Unfold< Norm<Expr, 2> >::apply(expr);
    }

    // template<class Expr>
    // inline typename Fold< Norm<Expr, 2> >::Type shallow_fold(const Norm<Expr, 2> &expr)
    // {
    //     return Fold< Norm<Expr, 2> >::apply(expr);
    // }

    template<class Expr, int Type>
    class Traits<Norm<Expr, Type> > : public Traits<Expr> {
    };

    /**
     * @ingroup reductions
     * @brief   \f$ || \cdot ||_{2} \f$
     */
    template<class Derived>
    inline Norm<Derived, 2> norm2(const Expression<Derived> &expr) {

        return Norm<Derived, 2>(expr.derived());
    }

    /**
     * @ingroup reductions
     * @brief   \f$ || \cdot ||_{1} \f$
     */
    template<class Derived>
    inline Norm<Derived, 1> norm1(const Expression<Derived> &expr) {

        return Norm<Derived, 1>(expr.derived());
    }

    /**
     * @ingroup reductions
     * @brief   \f$ || \cdot ||_{\infty} \f$
     */
    template<class Derived>
    inline Norm<Derived, INFINITY_NORM_TAG> norm_infty(const Expression<Derived> &expr) {

        return Norm<Derived, INFINITY_NORM_TAG>(expr.derived());
    }

    template<class Expr, int Type>
    inline Size size(const Norm<Expr, Type> &/*expr*/)
    {
        Size s(1);
        s.set(0, 1);
        return s;
    }

    //Derived types
    template<class Left, class Right, int NormType>
    using Distance = utopia::Norm<Binary<Left, Right, Minus>, NormType>;

}
#endif //UTOPIA_UTOPIA_NORM_HPP
