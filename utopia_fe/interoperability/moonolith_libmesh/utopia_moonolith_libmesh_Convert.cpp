#include "utopia_moonolith_libmesh_Convert.hpp"

namespace utopia {

    void convert_matrix(const moonolith::SparseMatrix<double> &in, USparseMatrix &out) {
        USizeType nnz = in.local_max_entries_x_col();

        USizeType n_local_rows = in.local_rows();
        USizeType n_local_cols = in.local_cols();

        // out = local_sparse(n_local_rows, n_local_cols, nnz);

        USizeType rows = in.rows();
        USizeType cols = in.cols();

        auto &&comm = out.comm();

        USizeType begin = 0;
        comm.exscan_sum(&n_local_rows, &begin, 1);

        UIndexArray d_nnz(n_local_rows, 0), o_nnz(n_local_rows, 0);

        assert(comm.rank() != 0 || begin == 0);

        Range r(begin, begin + n_local_rows);

        for (auto it = in.iter(); it; ++it) {
            const UScalar row = it.row() - begin;
            if (r.inside(it.col())) {
                ++d_nnz[row];
            } else {
                ++o_nnz[row];
            }
        }

        auto lo = layout(comm, n_local_rows, n_local_cols, rows, cols);
        out.sparse(lo, d_nnz, o_nnz);

        {
            Write<USparseMatrix> write(out);
            for (auto it = in.iter(); it; ++it) {
                out.set(it.row(), it.col(), *it);
            }
        }
    }

    template <int Dim>
    void ConvertMesh<libMesh::MeshBase, moonolith::Mesh<double, Dim>>::apply(const libMesh::MeshBase &in,
                                                                             moonolith::Mesh<double, Dim> &out) {
        out.clear();
        out.reserve(in.n_active_local_elem(), in.n_local_nodes());

        out.set_manifold_dim(in.mesh_dimension());

        std::unordered_map<libMesh::dof_id_type, moonolith::Integer> mapping;
        moonolith::Integer node_idx = 0;

        for (auto it = in.local_nodes_begin(); it != in.local_nodes_end(); ++it, ++node_idx) {
            const auto &node = **it;
            mapping[node.id()] = node_idx;
            make(node, out.add_node());
        }

        for (auto it = elements_begin(in); it != elements_end(in); ++it) {
            const auto &elem = **it;
            auto &e = out.add_elem();
            e.type = convert(elem.type());
            e.block = elem.subdomain_id();
            e.is_affine = elem.has_affine_map();
            e.global_idx = elem.unique_id();

            const std::size_t nn = elem.n_nodes();

            e.nodes.resize(nn);

            for (std::size_t i = 0; i < nn; ++i) {
                const auto node_id = elem.node_id(i);

                auto it = mapping.find(node_id);

                if (it == mapping.end()) {
                    e.nodes[i] = node_idx;
                    mapping[node_id] = node_idx++;
                    make(elem.node_ref(i), out.add_node());
                } else {
                    e.nodes[i] = it->second;
                }
            }
        }

        out.finalize();
    }

    template class ConvertMesh<libMesh::MeshBase, moonolith::Mesh<double, 1>>;
    template class ConvertMesh<libMesh::MeshBase, moonolith::Mesh<double, 2>>;
    template class ConvertMesh<libMesh::MeshBase, moonolith::Mesh<double, 3>>;

    //////////////////////////////////////////////////////////////////////////////////////////

    template <int Dim>
    void ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, Dim>>>::apply(
        const LibMeshFunctionSpace &in_space,
        moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &out) {
        int var_num = in_space.subspace_id();
        auto &in = in_space.mesh();
        auto &dof_map = in_space.dof_map();
        apply(in, dof_map, var_num, out);
    }

    template <int Dim>
    void ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, Dim>>>::apply(
        const libMesh::MeshBase &in,
        const libMesh::DofMap &dof_map,
        unsigned int var_num,
        moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &out) {
        convert(in, out.mesh());

        const auto n_local_dofs = dof_map.n_local_dofs();

        auto &out_dof_map = out.dof_map();
        out_dof_map.set_n_local_dofs(n_local_dofs);
        out_dof_map.set_n_dofs(dof_map.n_dofs());
        out_dof_map.set_max_nnz(max_nnz_x_row(dof_map));

        unsigned int sys_num = dof_map.sys_number();
        std::vector<libMesh::dof_id_type> dof_indices;

        const long n_elem = in.n_active_local_elem();
        out_dof_map.resize(n_elem);

        long idx = 0;
        long n_local_elems = in.n_active_local_elem();

        moonolith::Communicator comm(in.comm().get());
        comm.exscan(&n_local_elems, &idx, 1, moonolith::MPISum());

        for (long i = 0; i < n_elem; ++i) {
            out_dof_map.dof_object(i).element_dof = idx++;
        }

        auto fe_type = dof_map.variable(var_num).type();

        SizeType local_el_idx = -1;
        auto ebegin = elements_begin(in);
        auto eend = elements_end(in);
        for (auto it = ebegin; it != ebegin; ++it) {
            auto elem_ptr = *it;

            ++local_el_idx;

            auto &dof_object = out_dof_map.dof_object(local_el_idx);
            dof_object.global_idx = elem_ptr->id();
            dof_object.block = elem_ptr->subdomain_id();
            dof_object.type = convert(elem_ptr->type(), fe_type);

            if (fe_type.order == libMesh::CONSTANT) {
                dof_map.dof_indices(elem_ptr, dof_indices, var_num);
                auto nn = dof_indices.size();
                assert(nn == 1);

                dof_object.dofs.resize(nn);

                for (std::size_t i = 0; i < nn; ++i) {
                    dof_object.dofs[i] = dof_indices[i];
                }

            } else {
                const std::size_t nn = elem_ptr->n_nodes();
                dof_object.dofs.resize(nn);

                for (std::size_t i = 0; i < nn; ++i) {
                    const auto &node_ref = elem_ptr->node_ref(i);
                    const auto dof = node_ref.dof_number(sys_num, var_num, 0);
                    dof_object.dofs[i] = dof;
                }
            }
        }
    }

    template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 1>>>;
    template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 2>>>;
    template class ConvertFunctionSpace<LibMeshFunctionSpace, moonolith::FunctionSpace<moonolith::Mesh<double, 3>>>;
}  // namespace utopia
