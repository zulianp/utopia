#include "utopia_MeshParamSmoother.hpp"
#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"

namespace utopia {

    void MeshParamSmoother::read(Input &is)
    {
        is.get("operator-power", operator_power_);

        is.get("boundaries", [this](Input &in) {
            in.get_all([this](Input &in) {
                int id = -1;

                in.get("side", id);
                requested_boundary_ids_.insert(id);
            });
        });
    }

    static void constrain_p1_nodes_and_zero_out_p2(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        // LibMeshFunctionSpace &V,
        USparseMatrix &op,
        std::vector<UVector> &rhs)
    {
        // auto &mesh = V.mesh();
        // auto &dof_map = V.dof_map();
        const auto spatial_dim = mesh.spatial_dimension();

        auto n_local_dofs = dof_map.n_local_dofs();
        auto rr = row_range(op);

        UVector is_constrained = local_zeros(n_local_dofs);

        {
            Write<UVector> wi(is_constrained);

            for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
                auto &e = **e_it;

                for(unsigned int node_i = 0; node_i < e.n_nodes(); node_i++) {
                    const libMesh::Node * node = e.node_ptr(node_i);
                    const auto dof = node->dof_number(0, 0, 0);

                    auto idx = dof - rr.begin();

                    if(e.is_vertex(node_i)) {
                        is_constrained.set(dof, 1.);
                    }
                }
            }
        }

        for(auto &r : rhs) {
            r = e_mul(is_constrained, r);
        }

        set_zero_rows(op, is_constrained, 1.0);
    }

    static void assemble_laplacian(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, USparseMatrix &op)
    {
        //auto &mesh = V.mesh();
        //auto &dof_map = V.dof_map();
        const auto dim = mesh.mesh_dimension();

        op = sparse({dof_map.n_dofs(), dof_map.n_dofs()}, dof_map.get_n_nz(), dof_map.get_n_oz());

        auto fe = libMesh::FEBase::build(dim, dof_map.variable_type(0));

        int var_order = dof_map.variable_type(0).order;
        int quadrature_order = var_order + var_order;

        libMesh::QGauss q_gauss(dim, libMesh::Order(quadrature_order));
        fe->attach_quadrature_rule(&q_gauss);
        auto &dphi = fe->get_dphi();
        auto &JxW  = fe->get_JxW();

        libMesh::DenseMatrix<libMesh::Real> mat;
        std::vector<libMesh::dof_id_type> dof_indices;

        {
            Write<USparseMatrix> w_op(op, utopia::GLOBAL_ADD);

            for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
                auto &e = **e_it;
                dof_map.dof_indices(&e, dof_indices);

                auto n = dof_indices.size();

                mat.resize(n, n);
                mat.zero();
                fe->reinit(&e);

                auto nqp = JxW.size();
                auto n_shape_funcions = dphi.size();

                for(uint si = 0; si < n_shape_funcions; ++si) {
                    for(uint sj = 0; sj < n_shape_funcions; ++sj) {
                        for(uint k = 0; k < nqp; ++k) {
                            mat(si, sj) += (dphi[si][k] * dphi[sj][k]) * JxW[k];
                        }
                    }
                }

                add_matrix(mat, dof_indices, dof_indices, op);
            }
        }
    }

    static void assemble_laplacian(LibMeshFunctionSpace &V, USparseMatrix &op)
    {
        assemble_laplacian(V.mesh(), V.dof_map(), op);
    }

    static void fill_rhs(
        const libMesh::MeshBase &mesh,
        const SizeType n_local_dofs,
        std::vector<UVector> &rhs)
    {
        rhs.resize(n_local_dofs);

        for(auto &r : rhs) {
          r = local_zeros(n_local_dofs);
        }

        auto spatial_dim = mesh.spatial_dimension();

        Write<std::vector<UVector>> w(rhs);

        for(auto it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it) {
            auto * node  = *it;
            const auto dof = node->dof_number(0, 0, 0);

            for(int d = 0; d < spatial_dim; ++d) {
                rhs[d].set(dof, (*node)(d));
            }
        }
    }

    static void solve_and_compute_displacement(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        const USparseMatrix &lapl_k,
        const std::vector<UVector> &coords,
        const std::vector<UVector> &rhs,
        const USparseMatrix &permutation,
        std::vector<UVector> &disp)
    {
        const auto spatial_dim = mesh.spatial_dimension();

        Factorization<USparseMatrix, UVector> solver;
        solver.update(make_ref(lapl_k));

        std::vector<UVector> sol(spatial_dim);
        disp.resize(spatial_dim);
        for(int d = 0; d < spatial_dim; ++d) {
            sol[d] = rhs[d];
            solver.apply(rhs[d], sol[d]);
            disp[d] = permutation * (sol[d] - coords[d]);
        }

        // write("r.m", rhs[0]);
        // write("L.m", lapl_k);
        // write("P.m", permutation);
        // write("u.m", sol[0]);
    }

    void MeshParamSmoother::apply_aux(
        const libMesh::MeshBase &surf_mesh,
        const libMesh::DofMap &surf_dof_map,
        const USparseMatrix &permutation,
        const libMesh::DofMap &dof_map,
        libMesh::MeshBase &mesh
        )
    {
        const auto spatial_dim = surf_mesh.spatial_dimension();

        USparseMatrix op, lapl_k;
        assemble_laplacian(surf_mesh, surf_dof_map, op);

        auto n_local_dofs = surf_dof_map.n_local_dofs();

        std::vector<UVector> coords;
        fill_rhs(surf_mesh, n_local_dofs, coords);

        std::vector<UVector> rhs = coords;

        lapl_k = op;

        for(int i = 1; i < operator_power_; ++i) {
            lapl_k *= op;
        }

        // lapl_k *= -1.;

        constrain_p1_nodes_and_zero_out_p2(
               surf_mesh,
               surf_dof_map,
               lapl_k,
               rhs
        );

        std::vector<UVector> disp;
        //compute displacement
        solve_and_compute_displacement(
            surf_mesh,
            surf_dof_map,
            lapl_k,
            coords,
            rhs,
            permutation,
            disp);

        
        UVector gx = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());

        std::vector<PetscInt>    dofs(1);
        std::vector<PetscScalar> vals(1);

        for(int i = 0; i < spatial_dim; ++i) {
            gx = disp[i];
            synchronize(gx);

            Read<UVector> r(gx);
            for(auto it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it) {
                auto * node  = *it;
                dofs[0] = node->dof_number(0, 0, 0);

                gx.get(dofs, vals);
                (*node)(i) += vals[0];
            }
        }
    }

    void MeshParamSmoother::apply(libMesh::MeshBase &mesh)
    {
        const auto spatial_dim = mesh.spatial_dimension();
        const auto mesh_dim    = mesh.mesh_dimension();

        const bool is_surf = spatial_dim > mesh_dim;

        LibMeshFunctionSpace V(mesh, libMesh::LAGRANGE, libMesh::SECOND);
        V.initialize();

        LibMeshFunctionSpaceAdapter adapter;

        adapter.extract_surface_init(
            make_ref(mesh),
            V.dof_map(),
            0
        );

        auto &surf_dof_map = adapter.surf_dof_map();
        auto &surf_mesh = adapter.mesh();

        apply_aux(surf_mesh, surf_dof_map, *adapter.permutation(), V.dof_map(), mesh);
    }

}
