#include "utopia_LibMeshDofMapAdapterNew.hpp"
#include <cassert>

#define CHRONO_START()            \
    utopia::Chrono macro_chrono_; \
    macro_chrono_.start();
#define CHRONO_END(macro_name_)                                        \
    {                                                                  \
        macro_chrono_.stop();                                          \
        std::cout << macro_name_ << ":" << macro_chrono_ << std::endl; \
    }

namespace utopia {

    void ElementDofMapAdapter::make_permutation(const ElementDofMapAdapter &from,
                                                const ElementDofMapAdapter &to,
                                                USparseMatrix &mat) {
        CHRONO_START();

        using Scalar = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar> vals(1, 1.0);

        auto max_nnz = from.max_nnz_;
        assert(max_nnz > 0);
        mat = local_sparse(to.n_local_dofs_, from.n_local_dofs_, max_nnz);

        std::size_t n_elems = from.dof_map_.n_elements();

        assert(n_elems == to.dof_map_.n_elements());

        Write<USparseMatrix> w(mat, utopia::GLOBAL_INSERT);

        for (std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map_.dofs(e);
            const auto &to_dofs = to.dof_map_.dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to = to_dofs.size();

            assert(n_from == n_to);

            for (std::size_t i = 0; i < n_to; ++i) {
                irows[0] = to_dofs[i];
                icols[0] = from_dofs[i];
                mat.set_matrix(irows, icols, vals);
            }
        }

        CHRONO_END("ElementDofMapAdapter::make_permutation");
    }

    void ElementDofMapAdapter::make_vector_permutation(const int dim,
                                                       const ElementDofMapAdapter &from,
                                                       const ElementDofMapAdapter &to,
                                                       USparseMatrix &mat) {
        CHRONO_START();

        using Scalar = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar> vals(1, 1.0);

        auto max_nnz = from.max_nnz_;
        assert(max_nnz > 0);
        mat = local_sparse(to.n_local_dofs_, from.n_local_dofs_ * dim, max_nnz);

        std::size_t n_elems = from.dof_map_.n_elements();

        assert(n_elems == to.dof_map_.n_elements());

        Write<USparseMatrix> w(mat, utopia::GLOBAL_INSERT);

        for (std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map_.dofs(e);
            const auto &to_dofs = to.dof_map_.dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to = to_dofs.size();

            assert(n_from == n_to);

            for (std::size_t i = 0; i < n_to; ++i) {
                for (int d = 0; d < dim; ++d) {
                    irows[0] = to_dofs[i] + d;
                    icols[0] = from_dofs[i] * dim + d;
                    mat.set_matrix(irows, icols, vals);
                }
            }
        }

        CHRONO_END("ElementDofMapAdapter::make_vector_permutation");
    }

    void ElementDofMapAdapter::make_tensorize_permutation(const int dim,
                                                          const ElementDofMapAdapter &from,
                                                          const ElementDofMapAdapter &to,
                                                          USparseMatrix &mat) {
        using Scalar = Traits<USparseMatrix>::Scalar;
        using SizeType = Traits<USparseMatrix>::SizeType;

        std::vector<SizeType> irows(1), icols(1);
        std::vector<Scalar> vals(1, 1.0);

        auto max_nnz = from.max_nnz_;
        assert(max_nnz > 0);
        mat = local_sparse(to.n_local_dofs_ * dim, from.n_local_dofs_, max_nnz);

        std::size_t n_elems = from.dof_map_.n_elements();

        assert(n_elems == to.dof_map_.n_elements());

        Write<USparseMatrix> w(mat, utopia::GLOBAL_INSERT);

        for (std::size_t e = 0; e < n_elems; ++e) {
            const auto &from_dofs = from.dof_map_.dofs(e);
            const auto &to_dofs = to.dof_map_.dofs(e);

            const auto n_from = from_dofs.size();
            const auto n_to = to_dofs.size();

            assert(n_from == n_to);

            for (std::size_t i = 0; i < n_to; ++i) {
                for (int d = 0; d < dim; ++d) {
                    irows[0] = to_dofs[i] + d;
                    icols[0] = from_dofs[i];
                    mat.set_matrix(irows, icols, vals);
                }
            }
        }
    }

    void ElementDofMapAdapter::init(const libMesh::MeshBase &mesh, const libMesh::DofMap &dof_map, const int var_num) {
        comm_ = mesh.comm().get();
        n_local_dofs_ = dof_map.n_local_dofs();
        max_nnz_ = max_nnz_x_row(dof_map);

        // I do not think we need anything but 0 at the moment
        unsigned int comp = 0;
        unsigned int sys_num = dof_map.sys_number();
        std::vector<libMesh::dof_id_type> dof_indices;

        const auto n_elem = mesh.n_active_local_elem();
        dof_map_.resize(n_elem);

        SizeType local_el_idx = -1;

        for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
            ++local_el_idx;

            auto &dof_object = dof_map_.dof_object(local_el_idx);
            dof_object.global_idx = elem_ptr->id();
            dof_object.block = elem_ptr->subdomain_id();

            dof_map.dof_indices(elem_ptr, dof_indices, var_num);
            auto n_dofs_x_el = dof_indices.size();
            dof_object.dofs.resize(n_dofs_x_el);

            for (std::size_t i = 0; i < n_dofs_x_el; ++i) {
                dof_object.dofs[i] = dof_indices[i];
            }
        }
    }

    void ElementDofMapAdapter::init_surf_to_vol(libMesh::MeshBase &vol_mesh,
                                                const libMesh::DofMap &volume_dof_map,
                                                libMesh::MeshBase &surf_mesh,
                                                const libMesh::DofMap &surf_dof_map,
                                                const int var_num,
                                                const int surf_var_num) {
        CHRONO_START();

        comm_ = vol_mesh.comm().get();
        n_local_dofs_ = volume_dof_map.n_local_dofs();
        max_nnz_ = max_nnz_x_row(surf_dof_map);

        // I do not think we need anything but 0 at the moment
        unsigned int comp = 0;
        unsigned int surf_comp = 0;
        unsigned int sys_num = volume_dof_map.sys_number();
        unsigned int surf_sys_num = surf_dof_map.sys_number();

        std::vector<libMesh::dof_id_type> dof_surf;
        std::vector<libMesh::dof_id_type> dof_vol;

        const auto n_elem = surf_mesh.n_active_local_elem();
        dof_map_.resize(n_elem);

        SizeType local_el_idx = -1;

        for (const auto &b_elem : surf_mesh.active_local_element_ptr_range()) {
            ++local_el_idx;

            auto &dof_object = dof_map_.dof_object(local_el_idx);

            dof_object.global_idx = b_elem->id();
            // dof_object.element_dof = b_elem->id(); //FIXME not
            dof_object.block = b_elem->subdomain_id();

            const libMesh::Elem *v_elem = b_elem->interior_parent();

            surf_dof_map.dof_indices(b_elem, dof_surf);
            auto n_dofs_x_el = dof_surf.size();

            dof_object.dofs.resize(n_dofs_x_el);

            // loop through all nodes in each boundary element.
            for (unsigned int node = 0; node < b_elem->n_nodes(); node++) {
                // Node in boundary element.
                const libMesh::Node *b_node = b_elem->node_ptr(node);

                for (unsigned int node_id = 0; node_id < v_elem->n_nodes(); node_id++) {
                    // Nodes in interior_parent element.
                    const libMesh::Node *v_node = v_elem->node_ptr(node_id);

                    const auto v_dof = v_node->dof_number(sys_num, var_num, comp);

                    if (v_node->absolute_fuzzy_equals(*b_node, 1e-14)) {
                        // Global dof_index for node in BoundaryMesh
                        const auto b_dof = b_node->dof_number(surf_sys_num, surf_var_num, surf_comp);

                        auto dof_it = std::find(dof_surf.begin(), dof_surf.end(), b_dof);
                        assert(dof_it != dof_surf.end());

                        auto local_ind = std::distance(dof_surf.begin(), dof_it);
                        dof_object.dofs[local_ind] = v_dof;
                    }
                }
            }
        }

        long idx = 0;
        long n_local_elems = surf_mesh.n_active_local_elem();

        moonolith::Communicator comm(surf_mesh.comm().get());
        comm.exscan(&n_local_elems, &idx, 1, moonolith::MPISum());

        assert(idx < surf_mesh.n_active_elem());

        for (long i = 0; i < n_local_elems; ++i) {
            dof_map_.dof_object(i).element_dof = idx++;
        }

        CHRONO_END("ElementDofMapAdapter::init_surf_to_vol");
    }

    void ElementDofMapAdapter::make_element_node_map(ElementDofMapAdapter &out) const {
        CHRONO_START();

        const auto n_elems = dof_map_.n_elements();

        out.comm_ = comm_;
        out.dof_map_.resize(n_elems);
        out.n_local_dofs_ = 0;
        out.max_nnz_ = 0;

        for (std::size_t i = 0; i < n_elems; ++i) {
            SizeType n_e_dofs = dof_map_.dofs(i).size();
            out.max_nnz_ = std::max(out.max_nnz_, n_e_dofs);
            out.n_local_dofs_ += n_e_dofs;
        }

        SizeType dofs_begin = 0;
        out.comm_.exscan(&out.n_local_dofs_, &dofs_begin, 1, moonolith::MPISum());

        SizeType dof_id = dofs_begin;
        for (std::size_t i = 0; i < n_elems; ++i) {
            const auto &dof_in = dof_map_.dof_object(i);
            auto &dof_out = out.dof_map_.dof_object(i);
            dof_out = dof_in;

            for (auto &d : dof_out.dofs) {
                d = dof_id++;
            }
        }

        CHRONO_END("ElementDofMapAdapter::make_element_node_map");
    }

}  // namespace utopia
