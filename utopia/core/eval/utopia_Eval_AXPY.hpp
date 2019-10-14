#ifndef UTOPIA_UTOPIA_EVAL_AXPY_HPP
#define UTOPIA_UTOPIA_EVAL_AXPY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus>, Traits, Backend> {
    public:
        typedef utopia::Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus> Expr;
        typedef EXPR_TYPE(Traits, Expr) Result;

        inline static Result apply(const Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus> &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct(expr.right());
            result.axpy(expr.left().left(), Eval<Left, Traits>::apply(expr.left().right()));

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };


    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Binary<Left, Number<ScalarT>, Multiplies>, Right, Plus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Right)
        apply(const Binary< Binary<Left, Number<ScalarT>, Multiplies>, Right, Plus > &expr)
        {
            EXPR_TYPE(Traits, Right) result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct(Eval<Right, Traits>::apply(expr.right()));
            result.axpy(expr.left().right(), Eval<Left, Traits>::apply(expr.left().left()));

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    template<class LeftDerived, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Tensor<LeftDerived, Order>, Binary<Number<ScalarT>, Tensor<RightTensor, Order>, Multiplies>, Plus>, Traits, Backend> {
    public:
        typedef utopia::Tensor<LeftDerived, Order> Left;
        typedef utopia::Tensor<RightTensor, Order> Right;
        typedef utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus> Expr;

        inline static EXPR_TYPE(Traits, Left) apply(const Expr &expr)
        {
            EXPR_TYPE(Traits, Left) result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct( Eval<Left, Traits>::apply(expr.left()) );
            result.axpy(
                expr.right().left(),
                Eval<Right, Traits>::apply(expr.right().right())
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    template<class LeftDerived, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Tensor<LeftDerived, Order>, Binary<Tensor<RightTensor, Order>, Number<ScalarT>, Multiplies>, Plus>, Traits, Backend> {
    public:
        typedef utopia::Tensor<LeftDerived, Order> Left;
        typedef utopia::Tensor<RightTensor, Order> Right;
        typedef utopia::Binary<Left, Binary<Right, Number<ScalarT>, Multiplies>, Plus> Expr;

        inline static EXPR_TYPE(Traits, Left) apply(const Expr &expr)
        {
            EXPR_TYPE(Traits, Left) result;

            UTOPIA_TRACE_BEGIN(expr);

            result.construct( Eval<Left, Traits>::apply(expr.left()) );
            
            result.axpy(
                expr.right().right(),
                Eval<Right, Traits>::apply(expr.right().left())
            );

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    //TODO
    // template<class Left, class Traits, int Backend>
    // class Eval<Binary<Left, Factory<Identity, 2>, Plus>, Traits, Backend> {
    // public:
    //     inline static typename TypeAndFill<Traits, Left>::Type apply(const Binary<Left, Factory<Identity, 2>, Plus> &expr) {
    //         static_assert(Left::Order == 2, "can only be instantiated for 2nd order tensors");

    //         UTOPIA_TRACE_BEGIN(expr);

    //         typename TypeAndFill<Traits, Left>::Type result = Eval<Left, Traits>::apply(expr.left());
    //         UTOPIA_BACKEND(Traits).mat_diag_shift(result, 1.0);


    //         UTOPIA_TRACE_END(expr);
    //         return result;
    //     }
    // };

    template<class Left, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Left)
        apply(const Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto result = Eval<Left, Traits>::apply(expr.left());
            result.shift_diag(expr.right().left());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    /// axpy specialization l = x + alpha * y
    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<
            Assign<Left,
                   Binary<Left,
                          Binary<Number<ScalarT>,
                                 Right,
                                 Multiplies>,
                          Plus
                         >
                    >, Traits, Backend> {
        public:
            typedef utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus> RExpr;
            typedef utopia::Assign<Left, RExpr> Expr;



        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&alpha = expr.right().right().left();
            auto &&l = Eval<Left, Traits, Backend>::apply(expr.left());
            auto &&y = Eval<Left, Traits, Backend>::apply(expr.right().left());
            auto &&x = Eval<Right, Traits, Backend>::apply(expr.right().right().right());

            if(l.same_object(y)) {
                l.axpy(alpha, x);
            } else if(l.same_object(x)) {
                l.scale(alpha);
                l.axpy(1.0, y);
            } else {
                l.construct(x);
                l.scale(alpha);
                l.axpy(1.0, y);
            }
            
            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    //axpy specialization
    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<
            Assign<Left,
                   Binary<Left,
                          Binary<Number<ScalarT>,
                                 Right,
                                 Multiplies>,
                          Minus
                         >
                    >, Traits, Backend> {
        public:
            typedef utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Minus> RExpr;
            typedef utopia::Assign<Left, RExpr> Expr;
            using Scalar = typename Traits::Scalar;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&res = Eval<Left, Traits, Backend>::apply(expr.left());
            auto &&l = Eval<Left, Traits, Backend>::apply(expr.right().left());
            auto &&r = Eval<Right, Traits, Backend>::apply(expr.right().right().right());
            Scalar alpha = expr.right().right().left();

            //res = l -alpha * r;
            if(res.same_object(l)) {
                res.axpy(-alpha, r);
            } else if(res.same_object(r)) {
                res.scale(-alpha);
                res.axpy(1.0, l);
            } else {
                res.construct(l);
                res.axpy(-alpha, r);
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    //axpy specialization
    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<
            Assign<Left,
                   Binary<Binary<Number<ScalarT>,
                                 Right,
                                 Multiplies>,
                          Left,
                          Plus
                         >
                    >, Traits, Backend> {
        public:
            typedef utopia::Binary<Binary<Number<ScalarT>, Right, Multiplies>, Left, Plus> RExpr;
            typedef utopia::Assign<Left, RExpr> Expr;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&alpha = expr.right().left().left();

            if(&expr.left() == &expr.right().right()) {
                auto &&l = Eval<Left, Traits, Backend>::apply(expr.left());

                l.axpy(
                    alpha,
                    Eval<Right, Traits, Backend>::apply(expr.right().left().right())
                );

            } else {
                auto &&l  = Eval<Left, Traits, Backend>::apply(expr.left());
                auto &&ll = Eval<Left, Traits, Backend>::apply(expr.right().right());
                auto &&rr = Eval<Right, Traits, Backend>::apply(expr.right().left().right());

                if(&rr == &l) {
                    expr.left().assign(Eval<RExpr, Traits, Backend>::apply(expr.right()));
                } else {
                    l.assign(ll);
                    l.axpy(alpha, rr);
                }
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    /// l = alpha * r1 - r2;
    template<class Left, typename ScalarT, class RFirst, class RSecond, class Traits, int Backend>
    class Eval<
            Assign<Left, Binary<Binary<Number<ScalarT>, RFirst, Multiplies>, RSecond, Minus>>,
            Traits,
            Backend> {
    public:
        template<class Expr>
        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);
            auto &l   = Eval<Left, Traits>::apply(expr.left());
            auto &&r1 = Eval<RFirst, Traits>::apply(expr.right().left().right());
            auto &&r2 = Eval<RSecond, Traits>::apply(expr.right().right());
            const ScalarT alpha = expr.right().left().left();

            if(l.same_object(r2)) {
                l.scale(-1.0);
                l.axpy(alpha, r1);
            } else if(l.same_object(r1)) {
                l.scale(alpha);
                l.axpy(-1.0, r2);
            } else {
                l.assign(r1);
                l.scale(alpha);
                l.axpy(-1.0, r2);
            }

            UTOPIA_TRACE_END_SPECIALIZED(expr);
        }

    };


    //FIXME WHY IS THIS NEVER INSTANTIATED?
    // template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    // class Eval<Binary<Binary<Left, Number<ScalarT>, Multiplies>, Right, Minus>, Traits, Backend> {
    // public:
    //     inline static EXPR_TYPE(Traits, Right)
    //     apply(const Binary<Binary<Left, Number<ScalarT>, Multiplies>, Right, Minus > &expr)
    //     {
    //         EXPR_TYPE(Traits, Right) result;

    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).zaxpy(
    //                 -expr.left().right(),
    //                  Eval<Left,  Traits>::apply(expr.left().left()),
    //                  Eval<Right, Traits>::apply(expr.right()),
    //                  result
    //         );

    //         UTOPIA_TRACE_END(expr);
    //         return result;
    //     }
    // };
}

#endif //UTOPIA_UTOPIA_EVAL_AXPY_HPP
