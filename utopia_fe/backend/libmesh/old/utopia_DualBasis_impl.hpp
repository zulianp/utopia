#ifndef UTOPIA_DUAL_BASIS_IMPL_HPP
#define UTOPIA_DUAL_BASIS_IMPL_HPP

#include "utopia_DualBasis.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_Utils.hpp"

namespace utopia {

    template <class EDofMap>
    bool DualBasis::build_global_trafo(const libMesh::MeshBase &mesh,
                                       const SizeType n_local_dofs,
                                       const EDofMap &dof_map,
                                       const UVector &elem_to_transform,
                                       const double alpha,
                                       USparseMatrix &mat,
                                       const bool inverse) {
        using SizeType = UTOPIA_SIZE_TYPE(USparseMatrix);

        auto e_begin = elements_begin(mesh);
        auto e_end = elements_end(mesh);

        SizeType nnz = 0;
        if (!dof_map.empty()) {
            // estimate
            nnz = dof_map[0].dofs.size() * 10;
        }

        mat = local_sparse(n_local_dofs, n_local_dofs, nnz);

        Read<UVector> r(elem_to_transform);
        Write<USparseMatrix> w(mat);

        std::vector<bool> is_node_on_boundary(n_local_dofs, false);

        auto rr = row_range(mat);

        libMesh::DenseMatrix<libMesh::Real> local_trafo, inv_trafo, local_trafo_t;
        std::vector<libMesh::dof_id_type> dofs;

        if (e_begin != e_end) {
            assemble_local_trafo((*e_begin)->type(), alpha, local_trafo, inv_trafo);

            if (inverse) {
                inv_trafo.get_transpose(local_trafo_t);
            } else {
                local_trafo.get_transpose(local_trafo_t);
            }
        }

        SizeType idx = 0;
        for (auto it = e_begin; it != e_end; ++it, ++idx) {
            const auto *e = *it;

            if (elem_to_transform.get(dof_map[idx].element_dof) > 0.) {
                // dof_map.dof_indices(e, dofs);

                const auto &dofs = dof_map[idx].dofs;
                set_matrix(local_trafo_t, dofs, dofs, mat);

                for (auto d : dofs) {
                    if (rr.inside(d)) {
                        is_node_on_boundary[d - rr.begin()] = true;
                    }
                }
            }
        }

        for (SizeType i = 0; i < n_local_dofs; ++i) {
            if (is_node_on_boundary[i]) {
                auto dof_I = i + rr.begin();
                mat.set(dof_I, dof_I, 1.);
            }
        }

        return true;
    }

}  // namespace utopia

#endif  // UTOPIA_DUAL_BASIS_IMPL_HPP
