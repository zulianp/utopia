//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_AXPY_HPP
#define UTOPIA_UTOPIA_EVAL_AXPY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Right)
        apply(const Binary<Binary<Number<ScalarT>, Left, Multiplies>, Right, Plus> &expr) {
            EXPR_TYPE(Traits, Right) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(result, Eval<Right, Traits>::apply(expr.right()) );
            UTOPIA_BACKEND(Traits).axpy(
                    result,
                    expr.left().left(),
                    Eval<Left, Traits>::apply(expr.left().right())
                    );

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

            UTOPIA_BACKEND(Traits).assign(result, Eval<Right, Traits>::apply(expr.right()) );
            UTOPIA_BACKEND(Traits).axpy(
                    result,
                    expr.left().right(),
                    Eval<Left, Traits>::apply(expr.left().left())
                    );


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    template<class LeftTensor, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Wrapper<LeftTensor, Order>, Binary<Number<ScalarT>, Wrapper<RightTensor, Order>, Multiplies>, Plus>, Traits, Backend> {
    public:
        typedef utopia::Wrapper<LeftTensor, Order> Left;
        typedef utopia::Wrapper<RightTensor, Order> Right;
        typedef utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus> Expr;
       
        inline static EXPR_TYPE(Traits, Left) apply(const Expr &expr)
        {
            EXPR_TYPE(Traits, Left) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(result, Eval<Left, Traits>::apply(expr.left()) );
            UTOPIA_BACKEND(Traits).axpy(
                    result,
                    expr.right().left(),
                    Eval<Right, Traits>::apply(expr.right().right())
                    );


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    ///hack that only allows for Terminal symbol on the left-hand side since the compiler thinks
    ///that this is ambiguous with Binary<Left, Right, Op>
    template<class LeftTensor, int Order, class RightTensor, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Wrapper<LeftTensor, Order>, Binary<Wrapper<RightTensor, Order>, Number<ScalarT>, Multiplies>, Plus>, Traits, Backend> {
    public:
        typedef utopia::Wrapper<LeftTensor, Order> Left;
        typedef utopia::Wrapper<RightTensor, Order> Right;
        typedef utopia::Binary<Left, Binary<Right, Number<ScalarT>, Multiplies>, Plus> Expr;
       
        inline static EXPR_TYPE(Traits, Left) apply(const Expr &expr)
        {
            EXPR_TYPE(Traits, Left) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).assign(result, Eval<Left, Traits>::apply(expr.left()) );
            UTOPIA_BACKEND(Traits).axpy(
                    result,
                    expr.right().right(),
                    Eval<Right, Traits>::apply(expr.right().left())
                    );


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    // template<class Left, typename ScalarT, class Traits, int Backend>
    // class Eval<Assign<Left,
    //                   Binary< Number<ScalarT>,
    //                           Factory<Identity, 2>,
    //                           Multiplies>
    //                  >,
    //            Traits, Backend> {
    // public:
    //     inline static void apply(const Assign<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies> > &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).build(
    //                Eval<Left, Traits>::apply(expr.left()),
    //                size(expr.right().right()),
    //                expr.right().right().type()
    //         );

    //         UTOPIA_BACKEND(Traits).scal(
    //                 expr.right().left(),
    //                 Eval<Left, Traits>::apply(expr.left()),
    //                 Eval<Left, Traits>::apply(expr.left())
    //         );

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    template<class Left, class Traits, int Backend>
    class Eval<Binary<Left, Factory<Identity, 2>, Plus>, Traits, Backend> {
    public:
        inline static typename TypeAndFill<Traits, Left>::Type apply(const Binary<Left, Factory<Identity, 2>, Plus> &expr) {
            static_assert(Left::Order == 2, "can only be instantiated for 2nd order tensors");

            UTOPIA_TRACE_BEGIN(expr);

            typename TypeAndFill<Traits, Left>::Type result = Eval<Left, Traits>::apply(expr.left());
            UTOPIA_BACKEND(Traits).mat_diag_shift(result, 1.0);


            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Left)
        apply(const Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            EXPR_TYPE(Traits, Left) result = Eval<Left, Traits>::apply(expr.left());
            UTOPIA_BACKEND(Traits).mat_diag_shift(result, expr.right().left());

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    template<class Left, class Right, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Binary<Left, Number<ScalarT>, Multiplies>, Right, Minus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Right)
        apply(const Binary<Binary<Left, Number<ScalarT>, Multiplies>, Right, Minus > &expr)
        {
            EXPR_TYPE(Traits, Right) result;

            UTOPIA_TRACE_BEGIN(expr);

            UTOPIA_BACKEND(Traits).zaxpy(
                    -expr.left().right(),
                     Eval<Left,  Traits>::apply(expr.left().left()),
                     Eval<Right, Traits>::apply(expr.right()),
                     result
            );

            UTOPIA_TRACE_END(expr);
            return result;
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
                          Plus
                         >
                    >, Traits, Backend> {
        public:
            typedef utopia::Binary<Left, Binary<Number<ScalarT>, Right, Multiplies>, Plus> RExpr;
            typedef utopia::Assign<Left, RExpr> Expr;

        inline static bool apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            if(&expr.left() == &expr.right().left()) {
                auto &&l = Eval<Left, Traits, Backend>::apply(expr.left());
                UTOPIA_BACKEND(Traits).axpy(
                    l,
                    expr.right().right().left(),
                    Eval<Right, Traits, Backend>::apply(expr.right().right().right())
                );
            } else {
                auto &&l  = Eval<Left, Traits, Backend>::apply(expr.left());
                auto &&ll = Eval<Left, Traits, Backend>::apply(expr.right().left());
                auto &&rr = Eval<Right, Traits, Backend>::apply(expr.right().right().right());
                
                auto &&alpha = expr.right().right().left();

                if(&rr == &l) {
                    // assert(false);
                    expr.left() = Eval<RExpr, Traits, Backend>::apply(expr.right());
                } else {
                    l = ll;
                    UTOPIA_BACKEND(Traits).axpy(
                        l,
                        alpha,
                        rr
                    );
                }
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

            if(&expr.left() == &expr.right().right()) {
                auto &&l = Eval<Left, Traits, Backend>::apply(expr.left());
                UTOPIA_BACKEND(Traits).axpy(
                    l,
                    expr.right().left().left(),
                    Eval<Right, Traits, Backend>::apply(expr.right().left().right())
                );
            } else {
                auto &&l  = Eval<Left, Traits, Backend>::apply(expr.left());
                auto &&ll = Eval<Left, Traits, Backend>::apply(expr.right().right());
                auto &&rr = Eval<Right, Traits, Backend>::apply(expr.right().left().right());
               
                auto &&alpha = expr.right().left().left();

                if(&rr == &l) {
                    // assert(false);
                    expr.left() = Eval<RExpr, Traits, Backend>::apply(expr.right());
                } else {
                    l = ll;
                    UTOPIA_BACKEND(Traits).axpy(
                        l,
                        alpha,
                        rr
                    );
                }
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_AXPY_HPP
