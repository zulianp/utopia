#ifndef UTOPIA_FE_EVAL_LOCAL_HPP
#define UTOPIA_FE_EVAL_LOCAL_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_MultiTensor.hpp"
#include "utopia_FormTensor.hpp"
#include "utopia_FiniteElement.hpp"


namespace utopia {

    //FIXME (specialize and move to libmesh backend)
    template<class Space, class Traits, int Backend>
    class Eval<
        Gradient<
            TrialFunction<
                FiniteElement<Space>
                >
            >, Traits, Backend> {
    public:

        using Expr = Gradient<
            TrialFunction<
                FiniteElement<Space>
                >
            >;

        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        template<class T>
        inline static void apply(const Expr &expr, FormTensor<T, 1> &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &fe = *expr.expr().space_ptr();
            const auto &ctx = fe.ctx();
            const auto &space = fe.space();
            const SizeType id = space.subspace_id();
            const auto &g = ctx.grad(id);

            assert(!g.empty());

            const SizeType n = g.size();
            const SizeType n_qp = g[0].size();

            if(result.size() != n) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                if(result[i].size() != n_qp) {
                    result[i].resize(n_qp);
                }

                for(SizeType k = 0; k < n_qp; ++k) {
                    result[i][k] = g[i][k];
                }
            }

            UTOPIA_TRACE_END(expr);
        }   

    };

    //FIXME (specialize and move to libmesh backend)
    template<class Space, class Traits, int Backend>
    class Eval<
            TrialFunction<
                FiniteElement<Space>
            >, Traits, Backend> {
    public:

        using Expr = TrialFunction<
                FiniteElement<Space>
            >;

        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        template<class T>
        inline static void apply(const Expr &expr, FormTensor<T, 0> &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &fe = *expr.space_ptr();
            const auto &ctx = fe.ctx();
            const auto &space = fe.space();
            const SizeType id = space.subspace_id();
            // const auto &f = ctx.fun(id);

            assert(!ctx.trial().empty());
            const auto &f = ctx.trial()[id]->get_phi();

            assert(!f.empty());

            const SizeType n = f.size();
            const SizeType n_qp = f[0].size();

            if(result.size() != n) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                if(result[i].size() != n_qp) {
                    result[i].resize(n_qp);
                }

                for(SizeType k = 0; k < n_qp; ++k) {
                    result[i][k] = f[i][k];
                }
            }

            UTOPIA_TRACE_END(expr);
        }   

    };

    template<class FE>
    class Measure : public Expression<Measure<FE>> {
    public:

        Measure(FE &fe) : fe_(fe) {}

        inline FE &fe() { return fe_; }
        inline const FE &fe() const { return fe_; }
    private:
        FE &fe_;
    };

    template<class FE>
    class Traits<Measure<FE>> : public Traits<FE> {};

    template<class FE>
    inline Measure<FE> measure(FE &fe)
    {
        return Measure<FE>(fe);
    }

    //FIXME (specialize and move to libmesh backend)
    template<class Space, class Traits, int Backend>
    class Eval<Measure<FiniteElement<Space>>, Traits, Backend> {
    public:
        using Expr = utopia::Measure<FiniteElement<Space>>;

        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        template<class T>
        inline static void apply(const Expr &expr, MultiTensor<T, 0> &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &fe = expr.fe();
            const auto &ctx = fe.ctx();
            const auto &dx = ctx.dx();

            const SizeType n = dx.size();

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                result[i] = dx[i];
            }

            UTOPIA_TRACE_END(expr);
        }
    };


    template<int N, class MeasureT, class... Args>
    class Form {};


    template<class MeasureT, class Left, class Right>
    class Form<2, MeasureT, Left, Right> : public Expression< Form<2, MeasureT, Left, Right> > {
    public:

        Form(const MeasureT &measure, const Left &l, const Right &r)
        : measure_(measure), left_(l), right_(r)
        {}

        inline const Left &left() const { return left_; }
        inline const Right &right() const { return right_; }
        inline const MeasureT &measure() const { return measure_; }

    private:
        const MeasureT &measure_;
        UTOPIA_STORE_CONST(Left)  left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template<class MeasureT, class Expr>
    class Form<1, MeasureT, Expr> : public Expression< Form<1, MeasureT, Expr> > {
    public:

        Form(const MeasureT &measure, const Expr &expr)
        : measure_(measure), expr_(expr)
        {}

        inline const Expr &expr() const { return expr_; }
        inline const MeasureT &measure() const { return measure_; }

    private:
        const MeasureT &measure_;
        UTOPIA_STORE_CONST(Expr)  expr_;
    };

    template<class MeasureT, class Left, class Right>
    Form<2, MeasureT, Left, Right> form2(const Left &left, const Right &right, const MeasureT &measure)
    {
        return Form<2, MeasureT, Left, Right>(measure, left, right);
    }

    template<class MeasureT, class Expr>
    Form<1, MeasureT, Expr> form1(const Expr &expr, const MeasureT &measure)
    {
        return Form<1, MeasureT, Expr>(measure, expr);
    }

    template<class T, int Order, class Traits, int Backend>
    class Eval< FormTensor<T, Order>, Traits, Backend> {
    public:
        inline static const FormTensor<T, Order> & apply(const FormTensor<T, Order> &expr)
        {
            return expr;
        }

        inline static FormTensor<T, Order> apply(FormTensor<T, Order> &&expr)
        {
            return std::move(expr);
        }
    };

    template<class MeasureT, class InnerExpr, class Traits, int Backend>
    class Eval<Form<1, MeasureT, InnerExpr>, Traits, Backend> {
    public:

        using Expr = utopia::Form<1, MeasureT, InnerExpr>;

        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        template<class T>
        inline static void apply(const Expr &expr, Tensor<T, 1> &t_result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&inner_expr = Eval<InnerExpr, Traits>::apply(expr.expr());
            const auto &dx = expr.measure();

            const SizeType rows = inner_expr.size();
            const SizeType n_qp = dx.size();

            auto &result = t_result.derived();

            if(rows != result.size()) {
                result.resize(rows);
            }

            result.set(0.0);
            for(SizeType r = 0; r < rows; ++r) {
                for(SizeType k = 0; k < n_qp; ++k) {
                    result.add(r, inner_expr[r][k] * dx[k]);
                }
            }

            UTOPIA_TRACE_END(expr);
        }

    };


    template<class MeasureT, class Left, class Right, class Traits, int Backend>
    class Eval<Form<2, MeasureT, Left, Right>, Traits, Backend> {
    public:

        using Expr = utopia::Form<2, MeasureT, Left, Right>;

        using SizeType = typename Traits::SizeType;
        using Scalar   = typename Traits::Scalar;

        template<class T>
        inline static void apply(const Expr &expr, Tensor<T, 2> &t_result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            auto &&left  = Eval<Left, Traits>::apply(expr.left());
            auto &&right = Eval<Right, Traits>::apply(expr.right());
            const auto &dx = expr.measure();

            const SizeType rows = left.size();
            const SizeType cols = right.size();
            const SizeType n_qp = dx.size();

            auto &result = t_result.derived();

            if(rows != result.rows() || cols != result.cols()) {
                result.resize(rows, cols);
            }

            result.set(0.0);

            for(SizeType r = 0; r < rows; ++r) {
                for(SizeType c = 0; c < cols; ++c) {
                    for(SizeType k = 0; k < n_qp; ++k) {
                        result.add(r, c, inner(left[r][k], right[c][k]) * dx[k]);
                    }
                }
            }

            UTOPIA_TRACE_END(expr);
        }

    };

}


#endif //UTOPIA_FE_EVAL_LOCAL_HPP
