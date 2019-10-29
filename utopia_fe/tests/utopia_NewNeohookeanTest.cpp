#include "utopia_NewNeohookeanTest.hpp"

#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_Flow.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"

namespace utopia {


    template<class Vector, class Traits, int Backend>
    class Eval<Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
                    >
                >, Traits, Backend> 
    {
    public:

        using Expr = utopia::Interpolate<
                Vector,
                TrialFunction<
                    FiniteElement<ProductFunctionSpace<LibMeshFunctionSpace>, LIBMESH_TAG> 
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

        template<class T1>
        static void gather_interp_values(
            const FiniteElement<ProductFunctionSpace<T1>> &fe,
            const UVector &c,
            USerialVector &element_values
            )
        {
            using IndexArray = typename utopia::Traits<UVector>::IndexArray;

            const auto &tp_space = fe.space();
            const auto &ctx = fe.ctx();
            const auto &sub_0 = tp_space.subspace(0);
            const auto &mesh = sub_0.mesh();
            const auto &dof_map  = sub_0.dof_map();

            const auto &elem_ptr = mesh.elem(ctx.current_element());

            IndexArray prod_indices;
            std::vector<libMesh::dof_id_type> indices;

            tp_space.each([&](const int sub_index, const LibMeshFunctionSpace &space) {
                dof_map.dof_indices(elem_ptr, indices, space.subspace_id());
                prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
            });

            const std::size_t n_indices = prod_indices.size();
            element_values = zeros(n_indices);

            Write<USerialVector> w(element_values);
            Read<UVector> r(c);
            assert( c.has_ghosts() || c.comm().size() == 1 );
            c.get(prod_indices, element_values.entries());
        }

        template<class T1, class T2>
        inline static void eval_grad(
            const FiniteElement<ProductFunctionSpace<T1>> &fe,
            const UVector &c,
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
    

    void NewNeohookeanTest::run(Input &in)
    {
        std::cout << "[NewNeohookeanTest]" << std::endl;
        using Space = utopia::LibMeshFunctionSpace;
        using FE = utopia::FiniteElement<Space>;
        using TFE = utopia::FiniteElement<ProductFunctionSpace<Space>>;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<Space> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<Space, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);

        auto &V =space.space();

        //FIXME
        double lambda = 1.0, mu = 1.0, rescaling = 1.0;

        //FIXME use ghost
        UVector x = local_zeros(V[0].dof_map().n_local_dofs());

        //linear
        FormMatrixd gu, gu_t, alpha_gu; 
        MultiScalard dx;

        //non-linear
        MultiMatrixd guk, F, F_inv, F_inv_t, P;
        MultiScalard J, log_J, lambda_log_J, beta;
        FormMatrixd FG, FGF, stress_lin;
        FormScalard inner_Fg, PdotG;
        
        USparseMatrix A;
        UVector rhs;

        assemble(
            V,
            inner(grad(trial(V)), grad(trial(V))) * dX,
            A,
            [&](TFE &element, USerialMatrix &mat)
            {
                //Symbolic
                auto u = trial(element);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                //Numeric
                F = identity() + guk;
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);

                J = det(F);
                log_J = logn(J);
                lambda_log_J = (rescaling * lambda) * log_J;

                //////////////////////////////////////////
                //Stress linearization

                FG = F_inv_t * gu;
                FGF = FG * F_inv_t;

                alpha_gu = (rescaling * mu) * gu;

                beta = (rescaling * (lambda  * log_J - mu));
                inner_Fg = inner(F_inv_t, gu);
                stress_lin = alpha_gu - beta * FGF + lambda * (inner_Fg * F_inv_t);
                
                //////////////////////////////////////////
                //Form assembly


                mat = form2(stress_lin, gu, dx);
            }
        );


        assemble(
            V,
            grad(trial(V)) * dX,
            rhs,
            [&](TFE &element, USerialVector &vec)
            {
                //Symbolic
                auto u = trial(element);
                
                //Symbolic to Numeric
                gu = grad(u);
                guk = grad(interpolate(x, u));
                dx = measure(element);

                //Numeric
                F = identity() + guk;
                F_inv = inv(F);
                F_inv_t = transpose(F_inv);

                J = det(F);
                log_J = logn(J);
                lambda_log_J = (rescaling * lambda) * log_J;
                
                //////////////////////////////////////////
                //First Piola-Kirchoff tensor

                //axpy with scalar a
                P = F - F_inv_t;
                P *= (rescaling * mu);

                //axpy with multi-scalar a
                P += lambda_log_J * F_inv_t;
                PdotG = inner(P, gu);
                //////////////////////////////////////////
                //Form assembly
                vec = form1(PdotG, dx);
            }
        );

    }

}