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

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).zaxpy(
                    expr.left().left(),
                    Eval<Left,  Traits>::apply(expr.left().right()),
                    Eval<Right, Traits>::apply(expr.right()),
                    result);

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, typename ScalarT, class Traits, int Backend>
    class Eval<Assign<Left,
                      Binary< Number<ScalarT>,
                              Factory<Identity, 2>,
                              Multiplies>
                     >,
               Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies> > &expr) {
            UTOPIA_LOG_BEGIN(expr);

            UTOPIA_BACKEND(Traits).build(
                   Eval<Left, Traits>::apply(expr.left()),
                   size(expr.right().right()),
                   expr.right().right().type()
            );

            UTOPIA_BACKEND(Traits).scal(
                    expr.right().left(),
                    Eval<Left, Traits>::apply(expr.left()),
                    Eval<Left, Traits>::apply(expr.left())
            );

            //FIXME error handling

            UTOPIA_LOG_END(expr);
            return true;
        }
    };

    template<class Left, class Traits, int Backend>
    class Eval<Binary<Left, Factory<Identity, 2>, Plus>, Traits, Backend> {
    public:
        inline static typename TypeAndFill<Traits, Left>::Type apply(const Binary<Left, Factory<Identity, 2>, Plus> &expr) {
            static_assert(Left::Order == 2, "can only be instantiated for 2nd order tensors");

            UTOPIA_LOG_BEGIN(expr);

            typename TypeAndFill<Traits, Left>::Type result = Eval<Left, Traits>::apply(expr.left());
            const bool ok = UTOPIA_BACKEND(Traits).mat_diag_shift(result, 1.0);
            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };

    template<class Left, typename ScalarT, class Traits, int Backend>
    class Eval<Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus>, Traits, Backend> {
    public:
        inline static EXPR_TYPE(Traits, Left)
        apply(const Binary<Left, Binary<Number<ScalarT>, Factory<Identity, 2>, Multiplies>, Plus> &expr)
        {
            UTOPIA_LOG_BEGIN(expr);

            EXPR_TYPE(Traits, Left) result = Eval<Left, Traits>::apply(expr.left());
            const bool ok = UTOPIA_BACKEND(Traits).mat_diag_shift(result, expr.right().left());
            assert(ok);

            UTOPIA_LOG_END(expr);
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

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).zaxpy(
                    expr.left().right(),
                    Eval<Left,  Traits>::apply(expr.left().left()),
                    Eval<Right, Traits>::apply(expr.right()),
                    result
            );

            assert(ok);

            UTOPIA_LOG_END(expr);
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

            UTOPIA_LOG_BEGIN(expr);

            const bool ok = UTOPIA_BACKEND(Traits).zaxpy(
                    -expr.left().right(),
                     Eval<Left,  Traits>::apply(expr.left().left()),
                     Eval<Right, Traits>::apply(expr.right())
            );

            assert(ok);

            UTOPIA_LOG_END(expr);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_AXPY_HPP
