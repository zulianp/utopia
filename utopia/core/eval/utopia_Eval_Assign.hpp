#ifndef UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
#define UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Eval_Binary.hpp"

#include <type_traits>

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
            
            // expr.left().construct(
            //     Eval<Right, Traits>::apply(expr.right())
            // );

            Eval<Right, Traits, Backend>::apply(expr.right(), expr.left().derived());

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    template<class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Unary<Right, Minus>>, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Left, Unary<Right, Minus>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = Eval<Left, Traits>::apply(expr.left());
            
            l.construct(
                Eval<Right, Traits>::apply(expr.right().expr())
            );

            l.scale(-1.0);

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


    template<class Left, typename T, class Right, class Op, class Traits, int Backend>
    class Eval<
            Assign<Left, Binary<Number<T>, Unary<Right, Op>, Multiplies>>,
            Traits,
            Backend
            > {
    public:
        using Scalar = typename Traits::Scalar;

        inline static void apply(const Assign<Left, Binary<Number<T>, Unary<Right, Op>, Multiplies>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            auto &l = Eval<Left, Traits>::apply(expr.left());
            
            l.construct(
                Eval<Right, Traits>::apply(expr.right().right().expr())
            );

            l.transform(expr.right().right().operation());

            const Scalar alpha = expr.right().left();
            l.scale(alpha);

            UTOPIA_TRACE_END(expr);
        }
    };

    //Assign<V, Unary<Diag<M>, Op>>
    template<class V, class M, class Op, class Traits, int Backend>
    class Eval< Assign<Tensor<V, 1>, Unary<Diag<M>, Op>>, Traits, Backend> {
    public:
        static void apply(const Assign<Tensor<V, 1>, Unary<Diag<M>, Op>> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            auto &v = Eval<Tensor<V, 1>, Traits>::apply(expr.left());
            auto &&m = Eval<M, Traits>::apply(expr.right().expr().expr());
            auto &&op = expr.right().operation();

            m.build_diag(v);
            v.transform(op);

            UTOPIA_TRACE_END(expr);
        }
    };



    template<class Matrix, class Vector>
    using MatVecMult = Multiply<Tensor<Matrix, 2>, Tensor<Vector, 1> >;

    template<class Matrix, class Vector, class Op>
    using MatVecOpVec = 
        utopia::Binary<MatVecMult<Matrix, Vector>, 
                       Tensor<Vector, 1>, 
                       Op>;

    template<class Matrix, class Vector, class Op>
    using AssignMatVecOpVec = utopia::Assign<
                                        Tensor<Vector, 1>,
                                        MatVecOpVec<Matrix, Vector, Op>>;


    /// y = A * x op b
    template<class Matrix, class Vector, class Op, class Traits, int Backend>
    class Eval<AssignMatVecOpVec<Matrix, Vector, Op>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<Vector, 1>;
        using T2 = utopia::Tensor<Matrix, 2>;
        using Expr = utopia::AssignMatVecOpVec<Matrix, Vector, Op>;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &y  = Eval<T1, Traits>::apply(expr.left());
            auto &&A = Eval<T2, Traits>::apply(expr.right().left().left());
            auto &&x = Eval<T1, Traits>::apply(expr.right().left().right());
            auto &&b = Eval<T1, Traits>::apply(expr.right().right());
            auto &&op = expr.right().operation();

            if(y.is_alias(x)) {
                //temporary here
                Vector temp = x;
                A.multiply(temp, y);

            } else {
                A.multiply(x, y);
            }

            EvalBinaryAux<T1>::apply(y, b, op, y);
            UTOPIA_TRACE_END_SPECIALIZED(expr);
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
            if(tl.is_alias(res)) {
                res.axpy(-1.0, tr);
                return;
            }

            if(tr.is_alias(res)) {
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
            if(tl.is_alias(res)) {
                res.axpy(1.0, tr);
                return;
            }

            if(tr.is_alias(res)) {
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

    template<class Vector, typename T, class Traits, int Backend>
    class Eval<Assign<Tensor<Vector, 1>, Unary<Tensor<Vector, 1>, Reciprocal<T>> >, Traits, Backend> {
    public:
        using T1   = utopia::Tensor<Vector, 1>;
        using Expr = utopia::Assign<Tensor<Vector, 1>, Unary<Tensor<Vector, 1>, Reciprocal<T>> >;

        inline static bool apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);
            
            Vector &l = expr.left().derived();
            const Vector &r = expr.right().expr().derived();

            if(!l.is_alias(r)) {
                l.construct(r);
            } 
            
            const T num = expr.right().operation().numerator();
            l.reciprocal(num);                              

            UTOPIA_TRACE_END_SPECIALIZED(expr);
            return true;
        }
    };   


    template<typename T, class Vector>
    using ScaledVecExpr = utopia::Binary<Number<T>, Tensor<Vector, 1>, Multiplies>;


    template<class Matrix, class Vector, typename T, class Op>
    using AssignScaledVecOpMatVecMult = 
    utopia::Assign<Tensor<Vector, 1>, utopia::Binary<ScaledVecExpr<T, Vector>, MatVecMult<Matrix, Vector>, Op>>;

    template<class Matrix, class Vector, class Op>
    class EvalAssignScaledVecOpMatVecMult {
    public:
        using T1 = utopia::Tensor<Vector, 1>;
        using T2 = utopia::Tensor<Matrix, 2>;
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;

        template<typename T>
        static void apply(
            const ScaledVecExpr<T, Vector> &left_expr,
            const Op &op,
            const Matrix &A,
            const Vector &x,
            Vector &result)
        {
            auto &&left = Eval<T1, Traits>::apply(left_expr.right());
            const Scalar alpha = left_expr.left();

            if(result.is_alias(left)) {
                Vector temp;
                A.multiply(x, temp);
                result.scale(alpha);
                EvalBinaryAux<T1>::apply(result, temp, op, result);

            } else if(result.is_alias(x)) {
                Vector temp;
                A.multiply(x, temp);

                EvalBinaryAux<T1>::apply(left, temp, op, result);
            } else {
                A.multiply(x, result);

                if(std::is_same<Op, Minus>::value) {
                    result.axpy(-alpha, left);
                    result.scale(-1.0);
                } else if(std::is_same<Op, Plus>::value) {
                    result.axpy(alpha, left);
                } else {
                    assert(false && "Should not come here at the moment");
                    result = Binary<ScaledVecExpr<T, Vector>, Tensor<Vector, 1>, Op>(left_expr, result);
                }
            }
        }
    };

    /// result = fun(left) op A * right
    template<class Matrix, class Vector, typename T, class Traits, int Backend>
    class Eval<AssignScaledVecOpMatVecMult<Matrix, Vector, T, Plus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<Vector, 1>;
        using T2 = utopia::Tensor<Matrix, 2>;
        using Scalar = typename Traits::Scalar;

        using Expr = utopia::AssignScaledVecOpMatVecMult<Matrix, Vector, T, Plus>;
        static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &result  = Eval<T1, Traits>::apply(expr.left());
            auto &&A = Eval<T2, Traits>::apply(expr.right().right().left());
            auto &&right = Eval<T1, Traits>::apply(expr.right().right().right());

            auto &&op = expr.right().operation();

            EvalAssignScaledVecOpMatVecMult<Matrix, Vector, Plus>::apply(
                expr.right().left(),
                op,
                A,
                right,
                result
            );

            UTOPIA_TRACE_END_SPECIALIZED(expr);
        }
    };


    /// result = fun(left) op A * right
    template<class Matrix, class Vector, typename T, class Traits, int Backend>
    class Eval<AssignScaledVecOpMatVecMult<Matrix, Vector, T, Minus>, Traits, Backend> {
    public:
        using T1 = utopia::Tensor<Vector, 1>;
        using T2 = utopia::Tensor<Matrix, 2>;
        using Scalar = typename Traits::Scalar;

        using Expr = utopia::AssignScaledVecOpMatVecMult<Matrix, Vector, T, Minus>;
        static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN_SPECIALIZED(expr);

            auto &result  = Eval<T1, Traits>::apply(expr.left());
            auto &&A = Eval<T2, Traits>::apply(expr.right().right().left());
            auto &&right = Eval<T1, Traits>::apply(expr.right().right().right());

            auto &&op = expr.right().operation();

            EvalAssignScaledVecOpMatVecMult<Matrix, Vector, Minus>::apply(
                expr.right().left(),
                op,
                A,
                right,
                result
            );

            UTOPIA_TRACE_END_SPECIALIZED(expr);
        }
    };


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

            if(l.is_alias(r)) {
                l.scale(alpha);
            } else {
                l.construct(r);
                l.scale(alpha);
            }

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };


        
    template<class Left, class RowPtr, class ColIndex, class Values, class Traits, int Backend>
    class Eval< Assign<
                        Left,
                        Factory< CRS<RowPtr, ColIndex, Values>, 2>>,
                        Traits,
                        Backend
                        > {
    public:
        inline static bool apply(const Assign<Left, Factory<CRS<RowPtr, ColIndex, Values>, 2>> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);
            auto && left = Eval<Left, Traits>::apply(expr.left());
            auto && crs  = expr.right().type();
            auto s       = expr.right().size();

            //FIXME
            left.crs_init(
                left.comm().get(),
                INVALID_INDEX,
                INVALID_INDEX,
                s.get(0),
                s.get(1),
                crs.rowPtr(),
                crs.cols(),
                crs.values()
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }
    };




}

#endif //UTOPIA_UTOPIA_EVAL_ASSIGN_HPP_HPP
