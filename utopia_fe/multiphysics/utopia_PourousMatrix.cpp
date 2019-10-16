#include "utopia_PourousMatrix.hpp"
#include "utopia_LibMeshDofMapAdapterNew.hpp"

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
            FEBackend<LIBMESH_TAG>::subtract(f[qp], g[qp], ret[qp]);
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
        in.get("max-refinements", max_refinements_);
    }


    static void make_permutation(
                const int dim,
                const ElementDofMapAdapter &from,
                const ElementDofMapAdapter &to,
                 USparseMatrix &mat)
    {
        using Scalar   = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar>   vals(1, 1.0);

        auto max_nnz = from.max_nnz(); assert(max_nnz > 0);
        mat = local_sparse(to.n_local_dofs(), from.n_local_dofs(), max_nnz);

        std::size_t n_elems = from.dof_map().n_elements();
        
        assert(n_elems == to.dof_map().n_elements());

        Write<USparseMatrix>  w(mat, utopia::GLOBAL_INSERT);

        for(std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map().dofs(e);
            const auto &to_dofs   = to.dof_map().dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to   = to_dofs.size();

            assert(n_from == n_to);

            for(std::size_t i = 0; i < n_to; ++i) {
                for(int d = 0; d < dim; ++d) {
                    irows[0] = to_dofs[i] + d;
                    icols[0] = from_dofs[i];
                    mat.set_matrix(irows, icols, vals);
                }
            }
        }
    }

    static void make_permutation(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map_from,
        const int &var_num_from, 
        const libMesh::DofMap &dof_map_to,
        const int &var_num_to, 
        USparseMatrix &mat,
        const int tensor_dim = 1)
    {
        ElementDofMapAdapter a_from, a_to;
          a_from.init(
            mesh,
            dof_map_from,
            var_num_from);

          a_to.init(
            mesh,
            dof_map_to,
            var_num_to);

        make_permutation(tensor_dim, a_from, a_to, mat);

        assert(size(mat).get(0) == dof_map_to.n_dofs());
        assert(size(mat).get(1) == dof_map_from.n_dofs());
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::init(FunctionSpaceT &V)
    {
        // grad_space_.clear();
        if(grad_space_.empty()) {
            int dim = V.mesh().spatial_dimension();
            auto &aux = V.equation_systems().add_system<libMesh::LinearImplicitSystem>("gradient-recover");
            system_num_ = aux.number();

            // if(aux.is_initialized()) {
            //     std::cout << "Aux " << system_num_ << "/" << V.equation_systems().n_systems() << std::endl;
            //     std::cout << "Aux_n_var = " << aux.n_vars() << std::endl;

            //     for(int i = 0; i < dim + 1; ++i) {
            //         grad_space_ *= LibMeshFunctionSpace(aux, i);
            //     }

            //     return;
            // }

            grad_space_ *= LibMeshFunctionSpace(aux, aux.add_variable("fun", libMesh::Order(V.order(0)), libMesh::LAGRANGE) );

            for(int i = 0; i < dim; ++i) {
                grad_space_ *= LibMeshFunctionSpace(aux, aux.add_variable("grad_" + std::to_string(i), libMesh::Order(V.order(0)), libMesh::LAGRANGE) );
            }

            error_var_num_ = aux.add_variable("error_estimate", libMesh::Order(0), libMesh::MONOMIAL);
            grad_space_ *= LibMeshFunctionSpace(aux, error_var_num_);
            grad_space_.subspace(0).initialize();
        } 




        std::cout << "Aux " << system_num_ << "/" << V.equation_systems().n_systems() << std::endl;
        std::cout << "Aux_n_var = " << V.equation_systems().get_system(system_num_).n_vars() << std::endl;   
    }

    template<class Matrix, class Vector>
    bool GradientRecovery<Matrix, Vector>::refine(const Mortar<Matrix, Vector> &mortar, FunctionSpaceT &V, const Vector &sol)
    {
        estimate_error(mortar, V, sol);
        return apply_refinement(V);
    }

    template<class Matrix, class Vector>
    bool GradientRecovery<Matrix, Vector>::apply_refinement(FunctionSpaceT &V)
   {
        if(n_refinements_ >= max_refinements_) return false;

        auto &m = V.mesh();


        Read<UVector> r(error_);
        auto el = elements_begin(m);
        auto el_end = elements_end(m);

        bool refined = false;
        for(; el != el_end; ++el) {
            auto &e = **el;
            auto idx = e.dof_number(system_num_, error_var_num_, 0);

            assert(idx < grad_space_.subspace(0).dof_map().n_dofs());
            assert(idx < size(error_).get(0));

            auto err = error_.get(idx);
            //max_local_error_*=n_refinements_;

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
        // refinement.test_level_one(true);
        refinement.clean_refinement_flags();

        V.equation_systems().reinit();

        n_refinements_ += refined;
        return refined;
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::append_error_estimate(FunctionSpaceT &V)
    {
        utopia::convert(all_values_, *grad_space_[0].equation_system().solution);
        grad_space_[0].equation_system().solution->close();
    }

    template<class Matrix, class Vector>
    void GradientRecovery<Matrix, Vector>::estimate_error(const Mortar<Matrix, Vector> &mortar, FunctionSpaceT &V, const Vector &sol)
    {
        Chrono chrono; chrono.start();

        init(V);

        USparseMatrix permutation, mortar_matrix;
        if(!mortar.empty()) {
            make_permutation(
                V.mesh(),
                V.dof_map(),
                0,
                grad_space_.subspace(0).dof_map(),
                0,
                permutation,
                grad_space_.n_subspaces() - 1
            );

            mortar_matrix = USparseMatrix(permutation * *mortar.mortar_matrix_without_slave_dofs() * transpose(permutation));

            auto s_mm = size(mortar_matrix);
            disp(s_mm);
        }

      
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
        UIndexArray ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);
        UVector p_interp_buff = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
        copy_values(V, sol, P, p_interp_buff);
        synchronize(p_interp_buff);

        assert(!has_nan_or_inf(p_interp_buff));

        auto p_interp = interpolate(p_interp_buff, p);
        auto grad_form = inner(grad(p_interp), v) * dX;

        auto mass_form = inner(u, v) * dX;

        UVector grad_ph, grad_p_projected, mass_vec, inv_mass_vec;
        USparseMatrix mass_mat;
        utopia::assemble(grad_form, grad_ph);
        utopia::assemble(mass_form, mass_mat);

        assert(!has_nan_or_inf(grad_ph));
        assert(!has_nan_or_inf(mass_mat));

        mass_vec = sum(mass_mat, 1);
        e_pseudo_inv(mass_vec, inv_mass_vec);

        if(!mortar.empty()) {
            grad_ph  = transpose(mortar_matrix) * grad_ph;
            // //FIXME
            mass_mat = USparseMatrix(transpose(mortar_matrix) * mass_mat * mortar_matrix);
            
            grad_p_projected = mortar_matrix * e_mul(grad_ph, inv_mass_vec);
        } else {
            grad_p_projected = e_mul(grad_ph, inv_mass_vec);
        }


        //UNCOMMENT ME once the bug is fixed
        Adaptivity a;
        USparseMatrix pre_constraint, post_constraint;
        auto &W_i = W.subspace(0);
        a.constraint_matrix(W_i, pre_constraint, post_constraint);
        grad_p_projected += post_constraint * grad_p_projected;

        assert(!has_nan_or_inf(grad_ph));
        assert(!has_nan_or_inf(mass_vec));
        assert(!has_nan_or_inf(inv_mass_vec));
        assert(!has_nan_or_inf(grad_p_projected));

        UVector grad_p_projected_buff = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
        grad_p_projected_buff = grad_p_projected;
        synchronize(grad_p_projected_buff);


        assert(!has_nan_or_inf(grad_p_projected_buff));

        auto grad_interp = interpolate(grad_p_projected_buff, u);
        auto error_form  = inner(norm2(grad_interp - grad(p_interp)), test_e) * dX;
        auto vol_form    = inner(coeff(1.), test_e) * dX;

        UVector vol;
        utopia::assemble(error_form, error_);
        utopia::assemble(vol_form,   vol);
        // error = e_div(error, vol);

        assert(!has_nan_or_inf(error_));

        all_values_ = p_interp_buff + grad_p_projected + error_;

        assert(!has_nan_or_inf(all_values_));


        chrono.stop();
        std::cout << chrono << std::endl;

        append_error_estimate(V);
    }   

    template class GradientRecovery<USparseMatrix, UVector>;

}