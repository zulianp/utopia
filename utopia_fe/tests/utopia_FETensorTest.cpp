#include "utopia_FETensorTest.hpp"

#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_Flow.hpp"

#include "utopia_MultiTensor.hpp"
#include "utopia_FormTensor.hpp"
#include "utopia_FiniteElement.hpp"
#include "utopia_libmesh_FiniteElement.hpp"


#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"


namespace utopia {
    template<class Traits, int Backend>
    class Eval<
        Gradient<
            TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
                >
            >, Traits, Backend> {
    public:

        using Expr = Gradient<
            TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
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


    template<class Traits, int Backend>
    class Eval<
            TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
            >, Traits, Backend> {
    public:

        using Expr = TrialFunction<
                FiniteElement<LibMeshFunctionSpace>
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

    template<class LocalAssembler, class InitExpr>
    void assemble(
        FunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        USparseMatrix &matrix,
        LocalAssembler assembler
        )
    {
        const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

        Chrono c;
        c.start();

        auto &V = V_w.derived();
        auto &mesh = V.mesh();
        auto &dof_map = V.dof_map();

       
        libMesh::DofConstraints constraints;

        if(!disable_adaptivity) {
            Adaptivity::compute_all_constraints(
                mesh,
                dof_map,
                constraints
            );
        }

        LibMeshAssembler::allocate_matrix(
            V.dof_map(),
            matrix
        );


        USerialMatrix mat;

        std::vector<libMesh::dof_id_type> dofs;
        FiniteElement<LibMeshFunctionSpace> element(V);

        {
            Write<USparseMatrix> w(matrix, utopia::GLOBAL_ADD);

            auto e_begin = mesh.active_local_elements_begin();

            //Which basis-function
            auto u = trial(element);
            for(auto e_it = e_begin; e_it != mesh.active_local_elements_end(); ++e_it) {
                element.set((*e_it)->id());

                if(e_it == e_begin) {
                    element.init(init_expr.derived());
                } else {
                    element.reinit(init_expr.derived());
                }

                assembler(element, mat);

                // if(element.ctx().has_assembled()) {
                    V.dof_map().dof_indices(*e_it, dofs);

                    if(!disable_adaptivity) {
                        Adaptivity::constrain_matrix(*e_it, dof_map, constraints, mat, dofs);
                    }


                    add_matrix(mat, dofs, dofs, matrix);
                // }
            }
        }

        c.stop();
        std::cout << "assemble NEW(" << mesh.n_active_local_elem() << "):\n" << c << std::endl;
    }

    template<class InitExpr, class LocalAssembler>
    void assemble(
        FunctionSpace<LibMeshFunctionSpace> &V_w,
        const Expression<InitExpr> &init_expr,
        UVector &vector,
        LocalAssembler assembler)
    {
        const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

        Chrono c;
        c.start();

        auto &V = V_w.derived();
        auto &mesh = V.mesh();
        auto &dof_map = V.dof_map();

        libMesh::DofConstraints constraints;

        if(!disable_adaptivity) {
            Adaptivity::compute_all_constraints(
                mesh,
                dof_map,
                constraints
            );
        }

        typename Traits<UVector>::IndexSet ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);

        UVector temp_vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);

        {
            Write<UVector> w(temp_vec, utopia::GLOBAL_ADD);

            USerialVector vec;

            std::vector<libMesh::dof_id_type> dofs;
            FiniteElement<LibMeshFunctionSpace> element(V);

            auto e_begin = mesh.active_local_elements_begin();

            //Which basis-function
            auto u = trial(element);
            for(auto e_it = e_begin; e_it != mesh.active_local_elements_end(); ++e_it) {
                element.set((*e_it)->id());

                if(e_it == e_begin) {
                    element.init(init_expr.derived());
                } else {
                    element.reinit(init_expr.derived());
                }

                assembler(element, vec);

                // if(element.ctx().has_assembled()) {
                    V.dof_map().dof_indices(*e_it, dofs);

                    if(!disable_adaptivity) {
                        Adaptivity::constrain_vector(*e_it, dof_map, constraints, vec, dofs);
                    }

                    add_vector(vec, dofs, temp_vec);
                // }
            }
        }

        if(Traits<UVector>::Backend == utopia::TRILINOS) {
            vector = 1. * temp_vec;
        } else {
            vector = std::move(temp_vec);
        }

        c.stop();
        std::cout << "assemble NEW(" << mesh.n_active_local_elem() << "):\n" << c << std::endl;
    }

    void FETensorTest::run(Input &in)
    {
        std::cout << "[FETensorTest]" << std::endl;
        using Space = utopia::LibMeshFunctionSpace;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<Space> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<Space, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);

        bool no_solve = false;
        in.get("no-solve", no_solve);

        auto &V = space.space().subspace(0);
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////

        //Global tensors
        USparseMatrix A;
        UVector rhs, x, forcing_term;

        //Local tensors
        FormVectord g;
        FormScalard f;
        FormVectord alpha_g;
        MultiScalard dx;
            
        utopia::assemble(
            V,
            //FIXME
            inner(grad(trial(V)), grad(test(V))) * dX,
            A,
            [&](FiniteElement<Space> &element, USerialMatrix &mat) {
                //Symbolic
                auto u = trial(element);

                dx  = measure(element);
                g   = grad(u);
                alpha_g = 0.5 * g;
                mat = form2(alpha_g, g, dx);
            }
        );

        utopia::assemble(
            V,
            //FIXME
            inner(coeff(0.0), test(V)) * dX,
            rhs,
            [&](FiniteElement<Space> &element, USerialVector &vec) {
                //Symbolic
                auto u = trial(element);

                //Symbolic to Numeric conversion
                f = u;

                f *= 0.0;
                dx  = measure(element);
                vec = form1(f, dx);
            }
        );

        if(no_solve) return;

        x = local_zeros(local_size(A).get(0));
        forcing_term = local_zeros(local_size(x));
        forcing_function.eval(x, forcing_term);
        rhs += forcing_term;

        apply_boundary_conditions(V, A, rhs);

        Factorization<USparseMatrix, UVector> fact;
        fact.describe(std::cout);
        fact.solve(A, rhs, x);

        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}
