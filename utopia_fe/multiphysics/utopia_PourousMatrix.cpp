#include "utopia_PourousMatrix.hpp"

namespace utopia {
    //FIXME move me to the backend
    inline static auto apply_gradient_diff(
        const Interpolate<UVector,TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>> &left,
        const GradInterpolate<UVector, TrialFunction<LibMeshFunctionSpace> > &right,
        // const Op &op,
        const Minus &op,
        AssemblyContext<LIBMESH_TAG> &ctx) -> decltype( FEBackend<LIBMESH_TAG>::grad(right.expr().coefficient(), right.expr().fun(), ctx) )
    {
        auto f = FEBackend<LIBMESH_TAG>::fun(left.coefficient(), left.fun(), ctx);
        auto g = FEBackend<LIBMESH_TAG>::grad(right.expr().coefficient(), right.expr().fun(), ctx);
        auto ret = g;

        std::size_t n_quad_points = g.size();

        for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
            FEBackend<LIBMESH_TAG>::subtract(f[qp], g[qp].implementation(), ret[qp].implementation());
        }

        return std::move(ret);
    }

    using GradientDiffExpr = 
    Binary<
    Interpolate<UVector, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace> > >,
    Gradient<Interpolate<UVector, TrialFunction<LibMeshFunctionSpace> > >,
    Minus
    >;

    template<class Traits, int Backend, int IsQuadData>
    class FEEval<GradientDiffExpr, Traits, Backend, IsQuadData> {
    public:

        inline static auto apply(const GradientDiffExpr &expr, AssemblyContext<Backend> &ctx) -> decltype( apply_gradient_diff(expr.left(), expr.right(), expr.operation(), ctx) )
        {
            return apply_gradient_diff(expr.left(), expr.right(), expr.operation(), ctx);
        }
    };

    template<class Matrix, class Vector>
    GradientRecovery<Matrix, Vector>::GradientRecovery()
    : system_num_(-1), error_var_num_(-1), max_local_error_(10), n_refinements_(0), max_refinements_(5)
    {}

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::read(Input &in)
    {
        in.get("max-local-error", max_local_error_);
        in.get("max_refinements", max_refinements_);
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::init(FunctionSpaceT &V)
    {
        if(grad_space_.empty()) {
            int dim = V.mesh().spatial_dimension();
            auto &aux = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("gradient-recover");
            system_num_ = aux.number();
            grad_space_ *= LibMeshFunctionSpace(aux, aux.add_variable("p", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );

            for(int i = 0; i < dim; ++i) {
                grad_space_ *= LibMeshFunctionSpace(aux, aux.add_variable("grad_" + std::to_string(i), libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
            }

            error_var_num_ = aux.add_variable("e", libMesh::Order(0), libMesh::MONOMIAL);
            grad_space_ *= LibMeshFunctionSpace(aux, error_var_num_);
        }

        grad_space_.subspace(0).initialize();
    }

    template<class Matrix, class Vector>
    bool GradientRecovery<Matrix, Vector>::refine(const Mortar<Matrix, Vector> &mortar, FunctionSpaceT &V, const Vector &sol)
    {
        if(n_refinements_ >= max_refinements_) return false;

        init(V);
        recover(mortar, V, sol);
        return apply_refinement(V, sol);
    }

    // template<class Matrix, class Vector>
    // bool GradientRecovery<Matrix, Vector>::refine(FunctionSpaceT &V, const Vector &sol)
    // {
    //     if(n_refinements_ >= max_refinements_) return false;

    //     init(V);
    //     recover(V, sol);
    //     return apply_refinement(V, sol);
    // }

    template<class Matrix, class Vector>
    bool GradientRecovery<Matrix, Vector>::apply_refinement(FunctionSpaceT &V, const Vector &sol)
   {
        auto &m = V.mesh();

        Read<UVector> r(error_);

        bool refined = false;
        for(auto e_it = elements_begin(m); e_it != elements_end(m); ++e_it) {
            auto &e = **e_it;
            auto idx = e.dof_number(system_num_, error_var_num_, 0);

            auto err = error_.get(idx);

            if(err > max_local_error_) {
                refined = true;
                e.set_refinement_flag(libMesh::Elem::REFINE);
            }
        }


        libMesh::MeshRefinement refinement(m);
        auto e_it = elements_begin(m);

        if(e_it != elements_end(m)) {
            (*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
        }

        refinement.make_flags_parallel_consistent();
        refinement.refine_elements();
        refinement.test_level_one(true);

        m.prepare_for_use();
        V.initialize();

        n_refinements_ += refined;

        append_error_estimate(V);
        return refined;
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::append_error_estimate(FunctionSpaceT &V)
    {
        utopia::convert(all_values_, *grad_space_[0].equation_system().solution);
        grad_space_[0].equation_system().solution->close();
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::recover(const Mortar<Matrix, Vector> &mortar, const FunctionSpaceT &V, const Vector &sol)
    {
        Chrono chrono; chrono.start();
        const int dim = V.mesh().spatial_dimension();

        auto &P = grad_space_.subspace(0);
        auto W = grad_space_.subspace(1, dim + 1);
        auto E = grad_space_.subspace(dim + 1); //error space

        auto p = trial(P);
        auto u = trial(W);
        auto v = test(W);
        auto trial_e = trial(E);
        auto test_e  = test(E);

        auto &dof_map = P.dof_map();

        UVector p_interp_buff = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
        copy_values(V, sol, P, p_interp_buff);
        synchronize(p_interp_buff);

        auto p_interp = interpolate(p_interp_buff, p);
        auto grad_form = inner(grad(p_interp), v) * dX;

        auto mass_form = inner(u, v) * dX;

        UVector grad_ph, mass_vec;
        USparseMatrix mass_mat;
        utopia::assemble(grad_form, grad_ph);
        utopia::assemble(mass_form, mass_mat);
        mass_vec = sum(mass_mat, 1);

        UVector grad_p_projected = e_mul(grad_ph, 1./mass_vec);
        UVector grad_p_projected_buff = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
        grad_p_projected_buff = grad_p_projected;
        synchronize(grad_p_projected_buff);

        auto grad_interp = interpolate(grad_p_projected_buff, u);
        auto error_form  = inner(norm2(grad_interp - grad(p_interp)), test_e) * dX;
        auto vol_form    = inner(coeff(1.), test_e) * dX;

        UVector vol;
        utopia::assemble(error_form, error_);
        std::cout << tree_format(error_form.getClass()) << std::endl;
        utopia::assemble(vol_form,   vol);
        // error = e_div(error, vol);

        all_values_ = p_interp_buff + grad_p_projected + error_;



        chrono.stop();
        std::cout << chrono << std::endl;
    }   

    template class GradientRecovery<USparseMatrix, UVector>;

}