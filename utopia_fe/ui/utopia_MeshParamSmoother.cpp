#include "utopia_MeshParamSmoother.hpp"

#include "utopia_ZeroRowsToIdentity.hpp"

namespace utopia {

    void MeshParamSmoother::read(Input &is)
    {
        is.get("operator-power", operator_power_);

    }

    void MeshParamSmoother::apply(libMesh::MeshBase &mesh)
    {
        const auto spatial_dim = mesh.spatial_dimension();
        const auto mesh_dim = mesh.mesh_dimension();

        const bool is_surf = spatial_dim > mesh_dim;

        LibMeshFunctionSpace V(mesh, libMesh::LAGRANGE, libMesh::SECOND);
        V.initialize();

        auto &dof_map = V.dof_map();

        std::vector<bool> node_is_bdry;
        std::vector<UVector> rhs(spatial_dim);

        const auto n_local_dofs = dof_map.n_local_dofs();

        for(auto &r : rhs) {
          r = local_zeros(n_local_dofs);
        }


        USparseMatrix op = sparse({dof_map.n_dofs(), dof_map.n_dofs()}, dof_map.get_n_nz(), dof_map.get_n_oz());

        auto rr = row_range(op);

        if(is_surf) {
            assert(false && "IMPLEMENT ME");
        } else {

            const auto dim = mesh.mesh_dimension();

            auto fe = libMesh::FEBase::build(dim, dof_map.variable_type(0));
            
            int var_order = dof_map.variable_type(0).order;
            int quadrature_order = var_order + var_order;

            libMesh::QGauss q_gauss(dim-1, libMesh::Order(quadrature_order));
            fe->attach_quadrature_rule(&q_gauss);
            auto &dphi = fe->get_dphi();
            auto &JxW  = fe->get_JxW();


            std::vector<bool> has_values(n_local_dofs, false);
            
            libMesh::DenseMatrix<libMesh::Real> mat;
            std::vector<libMesh::dof_id_type> dof_indices;

            {
                Write<USparseMatrix> w_op(op, utopia::GLOBAL_ADD);
                for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
                    auto &e = **e_it;

                    auto n_sides = e.n_sides();

                    dof_map.dof_indices(&e, dof_indices);
                        
                    auto n = dof_indices.size();

                    mat.resize(n, n);
                    mat.zero();

                    for(uint i = 0; i < n_sides; ++i) {
                        if(e.neighbor_ptr(i) == nullptr) {
                            // auto side_ptr = e.build_side_ptr(i);
                            //const int side_id = mesh.get_boundary_info().boundary_id(&e, i);
                            fe->reinit(&e, i);

                            auto nqp = JxW.size();
                            auto n_shape_funcions = dphi.size();

                            for(uint si = 0; si < n_shape_funcions; ++si) {
                                for(uint sj = 0; sj < n_shape_funcions; ++sj) {
                                    for(uint k = 0; k < nqp; ++k) {
                                        mat(si, sj) += (dphi[si][k] * dphi[sj][k]) * JxW[k];
                                    }
                                }
                            }
                        }
                    }

                    add_matrix(mat, dof_indices, dof_indices, op);

                    for(uint i = 0; i < n; ++i) {
                        for(uint j = 0; j < n; ++j) {
                            if(std::abs(mat(i, j)) > 1e-16) {
                                has_values[dof_indices[i] - rr.begin()] = true;
                                break;
                            }
                        }
                    }
                } 
            }

            zero_rows_to_identity(op, 0.);
            USparseMatrix lapl_k = op;

            for(int i = 1; i < operator_power_; ++i) {
                lapl_k *= op;
            }


            UVector indicator = local_zeros(n_local_dofs);

            {
                Write<UVector> wi(indicator);
                Write<std::vector<UVector>> wrhs(rhs);

                for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
                    auto &e = **e_it;

                    auto n_sides = e.n_sides();

                    for (unsigned int node_i = 0; node_i < e.n_nodes(); node_i++)
                    {
                        if(e.is_vertex(node_i)) {
                            const libMesh::Node * node = e.node_ptr(node_i);
                            const auto dof = node->dof_number(0, 0, 0);

                            auto idx = dof - rr.begin();
                            if(rr.inside(dof) && has_values[idx]) {
                                indicator.set(dof, 1.);

                                for(int d = 0; d < spatial_dim; ++d) {
                                    rhs[d].set(dof, (*node)(d));
                                }
                            }
                        }
                    }
                }
            }

            set_zero_rows(lapl_k, indicator, 1.);

            Factorization<USparseMatrix, UVector> solver;
            solver.update(make_ref(lapl_k));

            std::vector<UVector> sol(spatial_dim);

            UVector gx = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());

            for(int d = 0; d < spatial_dim; ++d) {
                solver.apply(rhs[d], sol[d]);

                gx = sol[d];

                std::vector<PetscInt>   dofs(1);
                std::vector<PetscScalar> vals(1);

                for(auto it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it) {
                    auto * node  = *it;
                    dofs[0] = node->dof_number(0, 0, 0);

                    gx.get(dofs, vals);
                    (*node)(d) = vals[0];
                }
            }
        }

        V.initialize();
        V.equation_systems().clear();
    }

}
