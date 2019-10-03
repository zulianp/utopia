#ifndef UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Tracer.hpp"

namespace utopia {

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, Order>, Tensor<Right, Order> >, Traits, Backend> {
    public:
        typedef utopia::Tensor<Left, Order> LeftExpr;
        typedef utopia::Tensor<Right, Order> RightExpr;
        typedef utopia::Assign<LeftExpr, RightExpr> Expr;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<LeftExpr, Traits>::apply(expr.left()).assign(
                Eval<RightExpr, Traits>::apply(expr.right())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Right>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Right> &expr) {
            UTOPIA_TRACE_BEGIN(expr);
            
            expr.left().construct(
                Eval<Right, Traits>::apply(expr.right())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    template<class Left, class L, class R, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Binary<Tensor<L, Order>, Tensor<R, Order>, Minus>>, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Left, Binary<Tensor<L, Order>, Tensor<R, Order>, Minus>>;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&l   = Eval<Left, Traits>::apply(expr.left());
            auto &&b_l = Eval<Tensor<L, Order>, Traits>::apply(expr.right().left());
            auto &&b_r = Eval<Tensor<R, Order>, Traits>::apply(expr.right().right());
            
            apply_aux(b_l, expr.right().operation(), b_r, l);

            UTOPIA_TRACE_END(expr);
            return true;
        }

        template<class TL, class TR, class Result>
        static void apply_aux(const TL &tl, const Minus &, const TR &tr, Result &res)
        {
            if(tl.same_object(res)) {
                res.axpy(-1.0, tr);
                return;
            }

            if(tr.same_object(res)) {
                res.scale(-1.0);
                res.axpy(1.0, tl);
                return;
            }

            res = tl;
            res.axpy(-1.0, tr);
        }
    };


    template<class Left, class L, class R, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Binary<Tensor<L, Order>, Tensor<R, Order>, Plus>>, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Left, Binary<Tensor<L, Order>, Tensor<R, Order>, Plus>>;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&l   = Eval<Left, Traits>::apply(expr.left());
            auto &&b_l = Eval<Tensor<L, Order>, Traits>::apply(expr.right().left());
            auto &&b_r = Eval<Tensor<R, Order>, Traits>::apply(expr.right().right());
            
            apply_aux(b_l, expr.right().operation(), b_r, l);

            UTOPIA_TRACE_END(expr);
            return true;
        }

        template<class TL, class TR, class Result>
        static void apply_aux(const TL &tl, const Plus &, const TR &tr, Result &res)
        {
            if(tl.same_object(res)) {
                res.axpy(1.0, tr);
                return;
            }

            if(tr.same_object(res)) {
                res.axpy(1.0, tl);
                return;
            }

            res = tl;
            res.axpy(1.0, tr);
        }

    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Unary<Right, Abs> >, Traits, Backend> {
    public:
        typedef utopia::Assign<Left, Unary<Right, Abs> > Expr;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left = Eval<Left,  Traits>::apply(expr.left());
            left.construct(
                Eval<Right, Traits>::apply( expr.right().expr() )
            );

            left.transform( expr.right().operation() );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    //saves vector allocations but it is slower ???
    // template<class Left, class Right, class Op, class Traits, int Backend>
    // class Eval<Assign<Left, Unary<Right, Op> >, Traits, Backend> {
    // public:
    //     typedef utopia::Assign<Left, Unary<Right, Op> > Expr;

    //     inline static bool apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);
    //         UTOPIA_BACKEND(Traits).apply_unary(Eval<Left,  Traits>::apply(expr.left()),
    //                                            expr.right().operation(),
    //                                            Eval<Right, Traits>::apply( expr.right().expr()) );

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };   

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign<Left, Transposed <Tensor<Right, 2> > >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Transposed <Tensor<Right, 2> > > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left  = Eval<Left,  Traits>::apply(expr.left());
            auto &&right = Eval<Tensor<Right, 2>, Traits>::apply(expr.right().expr());

            right.transpose(left);

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    template<class Tensor, typename T, class Traits, int Backend>
    class Eval< Assign<Tensor, Binary<Number<T>, Tensor, Multiplies> >, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Tensor, Binary<Number<T>, Tensor, Multiplies> >;
        using Scalar = typename Traits::Scalar;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&l = Eval<Tensor, Traits>::apply(expr.left());
            auto &&r = Eval<Tensor, Traits>::apply(expr.right().right());
            const Scalar alpha = expr.right().left();

            if(l.same_object(r)) {
                l.scale(alpha);
            } else {
                l.construct(r);
                l.scale(alpha);
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

}

#endif //UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
