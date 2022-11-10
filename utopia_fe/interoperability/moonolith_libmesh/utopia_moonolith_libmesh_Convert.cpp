#include "utopia_moonolith_libmesh_Convert.hpp"

#include "utopia_libmesh_RetroCompatibility.hpp"
#include "utopia_libmesh_TransferUtils.hpp"

#include "utopia_MaxRowNNZ.hpp"

#include <libmesh/dof_map.h>
#include <libmesh/mesh_base.h>

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

        auto ebegin = in.active_local_elements_begin();
        auto eend = in.active_local_elements_end();

        for (auto it = ebegin; it != eend; ++it) {
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
    void convert_libmesh_to_moonolith(const libMesh::MeshBase &in,
                                      const libMesh::DofMap &dof_map,
                                      unsigned int var_num,
                                      moonolith::FunctionSpace<moonolith::Mesh<double, Dim>> &out,
                                      unsigned int comp) {
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
        auto ebegin = in.active_local_elements_begin();
        auto eend = in.active_local_elements_end();

        for (auto it = ebegin; it != eend; ++it) {
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

    template void convert_libmesh_to_moonolith<1>(const libMesh::MeshBase &in,
                                                  const libMesh::DofMap &dof_map,
                                                  unsigned int var_num,
                                                  moonolith::FunctionSpace<moonolith::Mesh<double, 1>> &out,
                                                  unsigned int comp);

    template void convert_libmesh_to_moonolith<2>(const libMesh::MeshBase &in,
                                                  const libMesh::DofMap &dof_map,
                                                  unsigned int var_num,
                                                  moonolith::FunctionSpace<moonolith::Mesh<double, 2>> &out,
                                                  unsigned int comp);

    template void convert_libmesh_to_moonolith<3>(const libMesh::MeshBase &in,
                                                  const libMesh::DofMap &dof_map,
                                                  unsigned int var_num,
                                                  moonolith::FunctionSpace<moonolith::Mesh<double, 3>> &out,
                                                  unsigned int comp);

    template <int Dim>
    void extract_surface(const libMesh::MeshBase &in,
                         ::moonolith::Mesh<double, Dim> &out_mesh,
                         const std::vector<int> &tags) {
        const bool select_all = tags.empty();  // FOR them moment tags are ignored
        out_mesh.clear();

        std::unordered_map<libMesh::dof_id_type, ::moonolith::Integer> mapping;
        long n_local_elems = 0;
        long n_nodes = 0;

        for (const auto &elem_ptr : in.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                n_local_elems++;

                auto side_ptr = elem_ptr->build_side_ptr(i);

                const std::size_t n_side_nodes = side_ptr->n_nodes();

                for (std::size_t k = 0; k < n_side_nodes; ++k) {
                    auto node_id = side_ptr->node_id(k);

                    auto it = mapping.find(node_id);

                    if (it == mapping.end()) {
                        mapping[node_id] = n_nodes++;
                    }
                }
            }
        }

        out_mesh.resize(n_local_elems, n_nodes);

        for (const auto &m : mapping) {
            make(in.node_ref(m.first), out_mesh.node(m.second));
        }

        std::size_t elem_idx = 0;
        for (const auto &elem_ptr : in.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                ///////////////////////////// MESH ////////////////////////////////
                auto side_ptr = elem_ptr->build_side_ptr(i);
                auto &e = out_mesh.elem(elem_idx);

                e.type = convert(side_ptr->type());
                e.block = utopia::libmesh::boundary_id(in.get_boundary_info(), elem_ptr, i);
                e.is_affine = side_ptr->has_affine_map();
                // e.global_idx = side_ptr->id();

                const std::size_t n_side_nodes = side_ptr->n_nodes();

                e.nodes.resize(n_side_nodes);

                for (std::size_t k = 0; k < n_side_nodes; ++k) {
                    auto node_id = side_ptr->node_id(k);
                    auto it = mapping.find(node_id);
                    assert(it != mapping.end());
                    e.nodes[k] = it->second;
                }

                elem_idx++;
            }
        }

        out_mesh.set_manifold_dim(Dim - 1);
        out_mesh.finalize();
    }

    template void extract_surface<1>(const libMesh::MeshBase &in,
                                     ::moonolith::Mesh<double, 1> &out_mesh,
                                     const std::vector<int> &tags);

    template void extract_surface<2>(const libMesh::MeshBase &in,
                                     ::moonolith::Mesh<double, 2> &out_mesh,
                                     const std::vector<int> &tags);

    template void extract_surface<3>(const libMesh::MeshBase &in,
                                     ::moonolith::Mesh<double, 3> &out_mesh,
                                     const std::vector<int> &tags);

    template <int Dim>
    void extract_trace_space(const libMesh::MeshBase &in,
                             const libMesh::DofMap &dof_map,
                             unsigned int var_num,
                             ::moonolith::FunctionSpace<::moonolith::Mesh<double, Dim>> &out,
                             const std::vector<int> &tags,
                             unsigned int comp) {
        const auto n_local_dofs = dof_map.n_local_dofs();
        const bool select_all = tags.empty();  // FOR them moment tags are ignored

        auto &out_mesh = out.mesh();
        auto &out_dof_map = out.dof_map();

        out_dof_map.clear();
        out_mesh.clear();

        out_dof_map.set_n_local_dofs(n_local_dofs);
        out_dof_map.set_n_dofs(dof_map.n_dofs());
        out_dof_map.set_max_nnz(max_nnz_x_row(dof_map));

        std::unordered_map<libMesh::dof_id_type, ::moonolith::Integer> mapping;
        long n_local_elems = 0;
        long n_nodes = 0;

        for (const auto &elem_ptr : in.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                n_local_elems++;

                auto side_ptr = elem_ptr->build_side_ptr(i);

                const std::size_t n_side_nodes = side_ptr->n_nodes();

                for (std::size_t k = 0; k < n_side_nodes; ++k) {
                    auto node_id = side_ptr->node_id(k);

                    auto it = mapping.find(node_id);

                    if (it == mapping.end()) {
                        mapping[node_id] = n_nodes++;
                    }
                }
            }
        }

        out_mesh.resize(n_local_elems, n_nodes);
        out_dof_map.resize(n_local_elems);

        for (const auto &m : mapping) {
            make(in.node_ref(m.first), out_mesh.node(m.second));
        }

        auto fe_type = dof_map.variable(var_num).type();
        unsigned int sys_num = dof_map.sys_number();

        std::size_t elem_idx = 0;
        for (const auto &elem_ptr : in.active_local_element_ptr_range()) {
            const std::size_t n_sides = elem_ptr->n_sides();

            for (std::size_t i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                ///////////////////////////// MESH ////////////////////////////////
                auto side_ptr = elem_ptr->build_side_ptr(i);
                auto &e = out_mesh.elem(elem_idx);

                e.type = convert(side_ptr->type());
                e.block = utopia::libmesh::boundary_id(in.get_boundary_info(), elem_ptr, i);
                e.is_affine = side_ptr->has_affine_map();
                // e.global_idx = side_ptr->id();

                const std::size_t n_side_nodes = side_ptr->n_nodes();

                e.nodes.resize(n_side_nodes);

                for (std::size_t k = 0; k < n_side_nodes; ++k) {
                    auto node_id = side_ptr->node_id(k);
                    auto it = mapping.find(node_id);
                    assert(it != mapping.end());
                    e.nodes[k] = it->second;
                }

                ///////////////////////////// DOFMAP ///////////////////////////

                auto &dof_object = out_dof_map.dof_object(elem_idx);
                // dof_object.global_idx = side_ptr->id();
                dof_object.block = e.block;
                dof_object.type = convert(side_ptr->type(), fe_type);

                const std::size_t n_vol_nodes = elem_ptr->n_nodes();

                dof_object.dofs.clear();
                dof_object.dofs.reserve(n_nodes);

                for (std::size_t k = 0; k < n_side_nodes; ++k) {
                    const auto &surf_node = side_ptr->node_ref(k);
                    if (!surf_node.has_dofs()) break;

                    for (std::size_t j = 0; j < n_vol_nodes; ++j) {
                        const auto &vol_node = elem_ptr->node_ref(j);

                        if (vol_node.absolute_fuzzy_equals(surf_node, 1e-14)) {
                            const auto v_dof = vol_node.dof_number(sys_num, var_num, comp);

                            dof_object.dofs.push_back(v_dof);
                            break;
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////

                elem_idx++;
            }
        }

        out_mesh.set_manifold_dim(Dim - 1);
        out_mesh.finalize();

        ////////////////////////////////////////////////////////////////

        long idx = 0;
        ::moonolith::Communicator comm(in.comm().get());
        comm.exscan(&n_local_elems, &idx, 1, ::moonolith::MPISum());

        assert(idx < in.n_active_elem());

        for (long i = 0; i < n_local_elems; ++i) {
            out_dof_map.dof_object(i).element_dof = idx++;
        }
    }

    template void extract_trace_space<1>(const libMesh::MeshBase &in,
                                         const libMesh::DofMap &dof_map,
                                         unsigned int var_num,
                                         ::moonolith::FunctionSpace<::moonolith::Mesh<double, 1>> &out,
                                         const std::vector<int> &tags = std::vector<int>(),
                                         unsigned int comp);

    template void extract_trace_space<2>(const libMesh::MeshBase &in,
                                         const libMesh::DofMap &dof_map,
                                         unsigned int var_num,
                                         ::moonolith::FunctionSpace<::moonolith::Mesh<double, 2>> &out,
                                         const std::vector<int> &tags = std::vector<int>(),
                                         unsigned int comp);

    template void extract_trace_space<3>(const libMesh::MeshBase &in,
                                         const libMesh::DofMap &dof_map,
                                         unsigned int var_num,
                                         ::moonolith::FunctionSpace<::moonolith::Mesh<double, 3>> &out,
                                         const std::vector<int> &tags = std::vector<int>(),
                                         unsigned int comp);

}  // namespace utopia
