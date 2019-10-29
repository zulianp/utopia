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

        template<class T, int Order>
        inline static void apply(const Expr &expr, FormTensor<T, Order> &result)
        {
            UTOPIA_TRACE_BEGIN(expr);
            const auto &fe = *expr.expr().space_ptr();
            eval_grad(fe, result);
            UTOPIA_TRACE_END(expr);
        }

        template<class T1, class T2>
        inline static void eval_grad(const FiniteElement<T1> &fe, FormTensor<T2, 1> &result)
        {
            

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

            
        }   


        template<class T1, class T2>
        inline static void eval_grad(const FiniteElement<ProductFunctionSpace<T1>> &fe, FormTensor<T2, 2> &result)
        {
            const auto &ctx = fe.ctx();
            const auto &space = fe.space();
            const SizeType id = space.subspace(0).subspace_id();
            const auto &g = ctx.vector_fe()[id]->grad;

            assert(!g.empty());

            const SizeType n = g.size();
            const SizeType n_qp = g[0].size();

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                if(n_qp != result[i].size()) {
                    result[i].resize(n_qp);
                }

                for(SizeType k = 0; k < n_qp; ++k) {
                    result[i][k] = g[i][k];
                }
            }
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

        inline std::string get_class() const override { return "Measure"; }
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


        inline std::string get_class() const override { return "Form2<" + left().get_class() + ", " + right().get_class() + ">"; }

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

        inline std::string get_class() const override { return "Form1<" + expr().get_class() + ">"; }

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

    template<class Vector, class Space, int FEBackendTag, class Traits, int Backend>
    class Eval<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<Space>, FEBackendTag> 
                    >
                >, Traits, Backend> 
    {
    public:

        using Expr = utopia::Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<Space>, FEBackendTag> 
                    >
                >;

        template<class T, int Order>
        inline static void apply(const Expr &expr, MultiTensor<T, Order> &result)
        {
            const auto &fe = *expr.expr().space_ptr();
            eval_fun(fe, result);
        }

        template<class T1, class T2>
        inline static void eval_fun(const FiniteElement<T1> &fe, FormTensor<T2, 0> &result)
        {
            assert(false && "IMPLEMENT ME");
        }


        template<class T1, class T2>
        inline static void eval_fun(const FiniteElement<ProductFunctionSpace<T1>> &fe, FormTensor<T2, 1> &result)
        {
            assert(false && "IMPLEMENT ME");
        }

    };

    template<class Vector, class Space, class Traits, int Backend>
    class Eval<Gradient<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<Space> 
                    >
                >>, Traits, Backend> {
    public:
        using Expr = utopia::Gradient<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<Space> 
                    >
                >>;

        template<class Out>
        inline static void apply(const Expr &expr, Out &result)
        {
            const auto &fe = *expr.expr().fun().space_ptr();
            eval_grad(
                fe,
                expr.expr().coefficient(),
                result
            );
        }

        template<class T1, class VectorNonConst>
        static void gather_interp_values(
            const FiniteElement<ProductFunctionSpace<T1>> &fe,
            const VectorNonConst &c,
            USerialVector &element_values
            )
        {
            using IndexArray = typename utopia::Traits<Vector>::IndexArray;

            const auto &tp_space = fe.space();
            const auto &ctx = fe.ctx();
            const auto &sub_0 = tp_space.subspace(0);
            const auto &mesh = sub_0.mesh();
            const auto &dof_map  = sub_0.dof_map();

            const auto &elem_ptr = mesh.elem(ctx.current_element());

            IndexArray prod_indices;
            std::vector<libMesh::dof_id_type> indices;

            tp_space.each([&](const int sub_index, const T1 &space) {
                dof_map.dof_indices(elem_ptr, indices, space.subspace_id());
                prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
            });

            const std::size_t n_indices = prod_indices.size();
            element_values = zeros(n_indices);

            Write<USerialVector> w(element_values);
            Read<VectorNonConst> r(c);
            assert( c.has_ghosts() || c.comm().size() == 1 );
            c.get(prod_indices, element_values.entries());
        }

        template<class T1, class T2, class VectorNonConst>
        inline static void eval_grad(
            const FiniteElement<ProductFunctionSpace<T1>> &fe,
            const VectorNonConst &c,
            MultiTensor<T2, 2> &result)
        {
            const auto &ctx = fe.ctx();
            const auto &space = fe.space();
            const SizeType id = space[0].subspace_id();
            const auto &g = ctx.vector_fe()[id]->grad;

            USerialVector element_values;
            gather_interp_values(fe, c, element_values);

            const SizeType rows = space.n_subspaces();
            const SizeType cols = space.subspace(0).mesh().spatial_dimension();
            Size s{rows, cols};

            const SizeType n_shape_functions = g.size();
            const SizeType n_quad_points = g[0].size();

            if(n_quad_points != result.size()) {
                result.resize(n_quad_points);
            }

            for(SizeType i = 0; i < n_quad_points; ++i) {
                result[i].resize(rows, cols);
                result[i].set(0.0);
            }

            for(std::size_t i = 0; i < n_shape_functions; ++i) {
                for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
                    result[qp] += element_values.get(i) * g[i][qp];
                }
            }
        }
    };


    template<class Matrix, class Op, class Traits, int Backend>
    class Eval< Binary<SymbolicTensor<Identity, 2>, Matrix, Op>, Traits, Backend> {
    public:
        using Expr = utopia::Binary<SymbolicTensor<Identity, 2>, Matrix, Op>;

        template<class Out>
        inline static void apply(const Expr &expr, Out &result)
        {
            auto &&mat = Eval<Matrix, Traits>::apply(expr.right());

            if(!result.is_alias(mat)) {
                result = mat;
            }

            if(std::is_same<Op, Plus>::value) {
                result.shift_diag(1.0);
            } else if(std::is_same<Op, Minus>::value) {
                result.shift_diag(-1.0);
            } else {
                assert(false && "IMPLEMENT ME");
            }

        }
    };

    template<class T, int Order>
    class EvalBinaryAux<Tensor<MultiTensor<T, Order>, Order>> {
    public:
        using Scalar = typename Traits<T>::Scalar;

        using TEval = utopia::EvalBinaryAux<T>;

        static void apply(const MultiTensor<Scalar, 0> &left,  const MultiTensor<T, Order> &right, const Multiplies &, MultiTensor<T, Order>  &result)
        {
            if(result.is_alias(right)) {
                result.scale(left);
            } else {
                result = right;
                result.scale(left);
            }
        }

        static void apply(
            const Scalar &left,
            const MultiTensor<T, Order> &right,
            const Multiplies &,
            MultiTensor<T, Order> &result)
        {
            if(!result.is_alias(right)) {
                result = right;
            }

            result.scale(left);
        }


        template<class Matrix, class Vector>
        static void apply(
            const MultiTensor<Matrix, 2> &left,
            const MultiTensor<Vector, 1> &right,
            const Multiplies &,
            MultiTensor<Vector, 1>  &result)
        {
            if(result.is_alias(right)) {
                assert(false);
                MultiTensor<Vector, 1> temp = right;
                left.multiply(temp, result);
            } else {
                left.multiply(right, result);
            }
        }

        template<class Op>
        static void apply(
            const MultiTensor<T, Order> &left,
            const MultiTensor<T, Order> &right,
            const Op &op,
            MultiTensor<T, Order>  &result)
        {
            const SizeType n = left.size();
            assert(n == right.size());

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                TEval::apply(left[i], right[i], op, result[i]);
            }
        }
    };


    template<class T, int Order>
    class EvalBinaryAux<Tensor<FormTensor<T, Order>, Order>> {
    public:
        using Scalar = typename Traits<T>::Scalar;
        using MTEval = utopia::EvalBinaryAux< Tensor<MultiTensor<T, Order>, Order >>;

        template<class Left, class Right, class Op, class Result>
        static void apply(const Left &left, const Right &right, const Op &op, Result &result)
        {
            assert(false && "IMPLEMENT ME");
        }

        static void apply(
            const Scalar left,
            const FormTensor<T, Order> &right,
            const Multiplies &,
            FormTensor<T, Order> &result)
        {
            if(!result.is_alias(right)) {
                result = right;
            }

            result.scale(left);
        }

        template<class MTType, int MTOrder, class Op>
        static void apply(
            const MultiTensor<MTType, MTOrder> &left,
            const FormTensor<T, Order> &right,
            const Op &op,
            FormTensor<T, Order> &result)
        {
            assert(left.size() == right[0].size());

            const SizeType n = right.size();
            const SizeType n_qp = left.size();

            if(result.size() != n) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                MTEval::apply(left, right[i], op, result[i]);
            }
        }

        template<class MTType, int MTOrder, class Op>
        static void apply(
            const FormTensor<MTType, MTOrder> &left,
            const MultiTensor<T, Order> &right,
            const Op &op,
            FormTensor<T, Order> &result)
        {
            assert(right.size() == left[0].size());

            const SizeType n = left.size();
            const SizeType n_qp = right.size();

            if(result.size() != n) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                MTEval::apply(left[i], right, op, result[i]);
            }
        }

        static void apply(const FormTensor<Scalar, 0> &left,  const FormTensor<T, Order> &right, const Multiplies &, FormTensor<T, Order>  &result)
        {
            if(result.is_alias(right)) {
                result.scale(left);
            } else {
                result = right;
                result.scale(left);
            }
        }

        template<class Matrix>
        static void apply(
            const FormTensor<Matrix, 2> &left,
            const FormTensor<T, Order> &right,
            const Multiplies &,
            FormTensor<T, Order>  &result)
        {
            if(result.is_alias(right)) {
                assert(false);
                MultiTensor<T, Order> temp = right;
                left.multiply(temp, result);
            } else {
                left.multiply(right, result);
            }
        }

        template<class Op>
        static void apply(
            const FormTensor<T, Order> &left,
            const FormTensor<T, Order> &right,
            const Op &op,
            FormTensor<T, Order> &result)
        {
            const SizeType n = left.size();
            assert(n == right.size());

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                MTEval::apply(left[i], right[i], op, result[i]);
            }
        }
    };

    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval< Dot<
                    Tensor<MultiTensor<Left, Order>, Order>,
                    Tensor<FormTensor<Right, Order>, Order>
                >, Traits, Backend> {
    public:

        using Expr = utopia::Dot<
                    Tensor<MultiTensor<Left, Order>, Order>,
                    Tensor<FormTensor<Right, Order>, Order>
                >;

        using Scalar = typename Traits::Scalar;
        using Result = utopia::FormTensor<Scalar, 0>;

        static void apply(const Expr &expr, Result &result)
        {
            auto &&left = expr.expr().left().derived();
            auto &&right = expr.expr().right().derived();

            const SizeType n = right.size();

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                left.dot(right[i], result[i]);
            }
        }
    };
    

}


#endif //UTOPIA_FE_EVAL_LOCAL_HPP
