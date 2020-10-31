//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_DIAG_HPP
#define UTOPIA_UTOPIA_EVAL_DIAG_HPP

#include "utopia_Diag.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Multiply.hpp"

namespace utopia {

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Left, Diag<Right>>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Left, Diag<Right>>;
        using Result = EXPR_TYPE(Traits, Left);

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            static_assert(Right::Order == 1, "Right has to be a vector");
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left = Eval<Left, Traits>::apply(expr.left());

            if (!result.is_alias(left)) {
                result.construct(left);
            }

            aux_apply(Eval<Right, Traits>::apply(expr.right().expr()), result);

            UTOPIA_TRACE_END(expr);
        }

        template <class T, class DiagVector>
        inline static void aux_apply(const DiagVector &d_mat, Tensor<T, 1> &result) {
            result.derived().e_mul(d_mat);
        }

        template <class T, class DiagVector>
        inline static void aux_apply(const DiagVector &d_mat, Tensor<T, 2> &result) {
            result.derived().diag_scale_right(d_mat);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Multiply<Diag<Left>, Right>, Traits, Backend> {
    public:
        using Expr = utopia::Multiply<Diag<Left>, Right>;
        using Result = EXPR_TYPE(Traits, Right);

        inline static Result apply(const Expr &expr) {
            Result result;
            apply(expr, result);
            return result;
        }

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&r = Eval<Right, Traits>::apply(expr.right());
            auto &&l = Eval<Left, Traits>::apply(expr.left().expr());

            if (!result.is_alias(r)) {
                result.construct(r);
            }

            // assert(!result.is_alias(l));

            aux_apply(l, result);

            UTOPIA_TRACE_END(expr);
        }

        template <class T, class DiagVector>
        inline static void aux_apply(const DiagVector &d_mat, Tensor<T, 1> &result) {
            result.derived().e_mul(d_mat);
        }

        template <class T, class DiagVector>
        inline static void aux_apply(const DiagVector &d_mat, Tensor<T, 2> &result) {
            result.derived().diag_scale_left(d_mat);
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Assign /////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Diag<Right>>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Diag<Right>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).diag(Eval<Right, Traits>::apply(expr.right().expr()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 1>, Diag<Right>>, Traits, Backend> {
    public:
        using LeftTensor = utopia::Tensor<Left, 1>;
        inline static bool apply(const Assign<LeftTensor, Diag<Right>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Right, Traits>::apply(expr.right().expr()).build_diag(Eval<LeftTensor, Traits>::apply(expr.left()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 2>, Diag<Diag<Right>>>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Tensor<Left, 2>, Diag<Diag<Right>>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Tensor<Left, 2>, Traits>::apply(expr.left())
                .diag(Eval<Right, Traits>::apply(expr.right().expr().expr()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, 1>, Diag<Diag<Right>>>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Tensor<Left, 1>, Diag<Diag<Right>>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).assign(Eval<Left, Traits>::apply(expr.right().expr().expr()));

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };

    template <class T, class Traits, int Backend>
    class Eval<Diag<Tensor<T, 2>>, Traits, Backend> {
    public:
        typedef utopia::Tensor<T, 2> WT;
        using Expr = utopia::Diag<WT>;
        using Result = typename Traits::Vector;

        inline static Result apply(const Expr &expr) {
            Result result;

            UTOPIA_TRACE_BEGIN(expr);

            Eval<WT, Traits>::apply(expr.expr()).build_diag(result);

            UTOPIA_TRACE_END(expr);
            return result;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////// Specialized /////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

    template <class Left, class Right, class Traits, int Backend>
    class Eval<InPlace<Tensor<Left, 2>, Diag<Tensor<Right, 1>>, Plus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<Left, 2>, Diag<Tensor<Right, 1>>, Plus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = expr.left().derived();
            auto &r = expr.right().expr().derived();

            l.shift_diag(r);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<InPlace<Tensor<Left, 2>, Diag<Tensor<Right, 1>>, Minus>, Traits, Backend> {
    public:
        using Expr = utopia::InPlace<Tensor<Left, 2>, Diag<Tensor<Right, 1>>, Minus>;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = expr.left().derived();
            Right r = expr.right().expr().derived();
            r.scale(-1.0);
            l.shift_diag(r);

            UTOPIA_TRACE_END(expr);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_DIAG_HPP
