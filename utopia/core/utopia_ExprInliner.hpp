#ifndef UTOPIA_EXPR_INLINER_HPP
#define UTOPIA_EXPR_INLINER_HPP

#include "utopia_Range.hpp"
#include <cassert>

namespace utopia {
    template<class Expr>
    class ExprInliner {
    public:
        typedef utopia::Traits<Expr> Traits;
        typedef typename Traits::SizeType SizeType;
        typedef typename Traits::Scalar Scalar;

        inline static Scalar eval_at(const Factory<Identity, 2> &, const SizeType i, const SizeType j)
        {
            return i == j;
        }

        inline static Scalar eval_at(const SymbolicTensor<Identity, 2> &, const SizeType i, const SizeType j)
        {
            return i == j;
        }

        inline static Scalar eval_at(const Factory<LocalIdentity, 2> &, const SizeType i, const SizeType j)
        {
            return i == j;
        }

        template<int Order>
        inline static Scalar eval_at(const Factory<Zeros, Order> &, const SizeType, const SizeType /*j = 0*/)
        {
            return 0.0;
        }

        template<typename T>
        inline static Scalar eval_at(const Factory<Values<T>, 2> &expr, const SizeType, const SizeType /*j = 0*/)
        {
            return expr.type().value();
        }

        template<int Order>
        inline static Scalar eval_at(const Factory<Zeros, Order> &, const SizeType)
        {
            return 0.0;
        }

        template<typename T>
        inline static Scalar eval_at(const Factory<Values<T>, 1> &expr, const SizeType)
        {
            return expr.type().value();
        }

        template<class InnerExpr>
        inline static Scalar eval_at(const Negate<InnerExpr> &expr, const SizeType i, const SizeType j)
        {
            return -eval_at(expr.expr(), i, j);
        }

        template<class Left, class Right>
        inline static void eval(const Construct<Left, Right> &expr){
            eval(expr.right(), expr.left());
        }

        template<class Left, class Right>
        inline static void eval(const Construct<Number<Left>, Right> &expr){
            expr.left() = eval_at(expr.right(), 0, 0);
        }

        template<class Left, class Right>
        inline static void eval(const Construct<Number<Left>, Determinant<Right>> &expr){
            expr.left() = Eval<Determinant<Right>, Traits, Traits::Backend>::apply(expr);
        }

        template<class InnerExpr>
        inline static Scalar eval(const Determinant<InnerExpr> &expr) {
            return Eval<Determinant<InnerExpr>, Traits, Traits::Backend>::apply(expr);
        }

        template<class Left, class Right>
        inline static void eval(const Assign<Left, Right> &expr){
            eval(expr.right(), expr.left());
        }

        template<class Left, class Right>
        inline static void eval(const Assign<Left, Inverse<Right>> &expr){
            Eval<Inverse<Right>, Traits, Traits::Backend>::apply(expr.right(), expr.left());
        }

        template<class Left, class Right, class Operation>
        inline static void eval(const InPlace<Left, Right, Operation> &expr)
        {
            //FIXME connect to backend without tree transformation
            typedef utopia::Binary<Left, Right, Operation> TransformedExpr;

                   eval(
                    TransformedExpr(
                            expr.left(),
                            expr.right(),
                            expr.operation()
                    ),  expr.left());
        }

        ////////////////////////////////////////////////
        //DerivedTensor order 2
        ////////////////////////////////////////////////

        template<class Derived, class DerivedTensor>
        inline static void eval(const Expression<Derived> &expr_w, Tensor<DerivedTensor, 2> &result)
        {
            const Derived &expr = expr_w.derived();
            Size s = size(expr);
            result.derived().resize(s);

            Write<DerivedTensor> w(result.derived());

            Range rr = row_range(result);
            const SizeType cols = s.get(1);
            for(auto i = rr.begin(); i < rr.end(); ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    result.derived().set(i, j, eval_at(expr, i, j));
                }
            }
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Left, Right, Operation> &expr, const SizeType i, const SizeType j)
        {
            return Operation::template apply<Scalar>(eval_at(expr.left(),  i, j),
                eval_at(expr.right(), i, j) );
        }

        template<class DerivedTensor>
        inline static Scalar eval_at(const Tensor<DerivedTensor, 2> &expr, const SizeType i, const SizeType j)
        {
            return expr.derived().get(i, j);
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Number<Left>, Right, Operation> &expr, const SizeType i, const SizeType j)
        {
            return Operation::template apply<Scalar>(static_cast<Left>(expr.left()),
                eval_at(expr.right(), i, j) );
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Left, Number<Right>, Operation> &expr, const SizeType i, const SizeType j)
        {
            return Operation::template apply<Scalar>(eval_at(expr.left(), i, j),
                static_cast<Right>(expr.right()));
        }

        template<class InnerExpr, class Operation>
        inline static Scalar eval_at(const Unary<InnerExpr, Operation> &expr, const SizeType i, const SizeType j)
        {
            return Operation::template apply<Scalar>(eval_at(expr.expr(), i, j));
        }

        template<class InnerExpr>
        inline static Scalar eval_at(const Unary<InnerExpr, Minus> &expr, const SizeType i, const SizeType j)
        {
            return -eval_at(expr.expr(), i, j);
        }

        template<class InnerExpr>
        inline static Scalar eval_at(const Transposed<InnerExpr> &expr, const SizeType i, const SizeType j)
        {
            return eval_at(expr.expr(), j, i);
        }

        ////////////////////////////////////////////////
        //DerivedTensors order 1
        ////////////////////////////////////////////////

        template<class Derived, class DerivedTensor>
        inline static void eval(const Expression<Derived> &expr_w, Tensor<DerivedTensor, 1> &result)
        {
            const Derived &expr = expr_w.derived();
            auto &d = result.derived();

            auto s = size(expr);
            d.resize(s);

            Write<DerivedTensor> w(d);

            Range r = range(result);
            for(auto i = r.begin(); i < r.end(); ++i) {
                d.set(i, eval_at(expr, i));
            }
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Left, Right, Operation> &expr, const SizeType i)
        {
            return Operation::template apply<Scalar>(eval_at(expr.left(),  i),
                eval_at(expr.right(), i) );
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Number<Left>, Right, Operation> &expr, const SizeType i)
        {
            return Operation::template apply<Scalar>(expr.left(),
                eval_at(expr.right(), i) );
        }

        template<class Left, class Right, class Operation>
        inline static Scalar eval_at(const Binary<Left, Number<Right>, Operation> &expr, const SizeType i)
        {
            return Operation::template apply<Scalar>(eval_at(expr.left(),  i),
                expr.right());
        }

        template<class InnerExpr, class Operation>
        inline static Scalar eval_at(const Unary<InnerExpr, Operation> &expr, const SizeType i)
        {
            return Operation::template apply<Scalar>(eval_at(expr.expr(),  i));
        }

        template<class DerivedTensor>
        inline static Scalar eval_at(const Tensor<DerivedTensor, 1> &expr, const SizeType i, const SizeType j)
        {
            UTOPIA_UNUSED(j);
            assert(j == 0 && "Trying to access tensor of order 1 like an order 2 one.");
            return expr.derived().get(i);
        }

        template<class DerivedTensor>
        inline static Scalar eval_at(const Tensor<DerivedTensor, 1> &expr, const SizeType i)
        {
            return expr.derived().get(i);
        }

        template<class Left, class Right>
        inline static Scalar eval_at(const Multiply<Left, Right> &expr,
            const SizeType i, const SizeType j) {

            //Not implemented in parallel
            //FIXME implement range(expr)

            Scalar result = 0;
            Size s = size(expr.left());
            const SizeType cols = s.get(1);

            assert(size(expr).n_dims() == 2 || j == 0);

            for(SizeType k = 0; k < cols; ++k) {
                result += eval_at(expr.left(), i, k) * eval_at(expr.right(), k, j);
            }

            return result;
        }

        template<class Left, class Right>
        inline static Scalar eval_at(const Multiply<Left, Right> &expr,
            const SizeType i) {

            //Not implemented in parallel
            //FIXME implement range(expr)

            Scalar result = 0;
            Size s = size(expr.right());
            const SizeType rows = s.get(0);
            for(SizeType k = 0; k < rows; ++k) {
                result += eval_at(expr.left(), i, k) * eval_at(expr.right(), k);
            }

            return result;
        }

        ////////////////////////////////////////////////
        //DerivedTensors order 0
        ////////////////////////////////////////////////

        template<typename T>
        inline static void eval(const Expr &expr, Number<T> &result)
        {
            result = eval_at(expr, 0, 0);
        }

        template<class InnerExpr, class ReduceOp, typename T, class Operation>
        inline static Scalar eval_at(const Binary< Reduce<InnerExpr, ReduceOp>, Number<T>, Operation> &expr, const int i = 0, const int j = 0)
        {
            assert(i == 0 && j == 0);
            return Operation::template apply<Scalar>(static_cast<Scalar>(eval_at(expr.left())), static_cast<T>(expr.right()));
        }

        template<class InnerExpr, class ReduceOp, typename T, class Operation>
        inline static Scalar eval_at(const Binary<Number<T>, Reduce<InnerExpr, ReduceOp>, Operation> &expr, const int i = 0, const int j = 0)
        {
            assert(i == 0 && j == 0);
            return Operation::template apply<Scalar>(static_cast<T>(expr.left()), static_cast<Scalar>(eval_at(expr.right())));
        }

        template<class Left, class ReduceOpLeft, class Right, class ReduceOpRight, class Operation>
        inline static Scalar eval_at(const Binary<Reduce<Left, ReduceOpLeft>, Reduce<Right, ReduceOpRight>, Operation> &expr, const int i = 0, const int j = 0)
        {
            assert(i == 0 && j == 0);
            return Operation::template apply<Scalar>(eval_at(expr.left()), eval_at(expr.right()));
        }

        template<class InnerExpr, class Operation>
        inline static Scalar eval_at(const Reduce<InnerExpr, Operation> &expr, const int i = 0, const int j = 0)
        {
            assert(i == 0 && j == 0);

            Size s = size(expr.expr());

            if(s.get(0) == 0) {
                return 0;
            }

            const SizeType rows = s.get(0);
            const SizeType cols = s.n_dims() == 1? 1 : s.get(1);

            Scalar result = eval_at(expr.expr(), 0, 0);
            if(s.n_dims() == 1) {
                for(SizeType i = 1; i < rows; ++i) {
                    result = Operation::template apply<Scalar>(result, eval_at(expr.expr(), i, 0));
                }
            } else {

                SizeType j = 1;
                for(SizeType i = 0; i < rows; ++i) {
                    for(; j < cols; ++j) {
                     result = Operation::template apply<Scalar>(result, eval_at(expr.expr(), i, j));
                    }

                    j = 0;
                }
            }

            return result;
        }

    };

    template<class Expr, class Result>
    void inline_eval(const Expression<Expr> &expr, Result &result)
    {
        ExprInliner<Expr>::eval(expr.derived(), result);
    }
}

#endif //UTOPIA_EXPR_INLINER_HPP
