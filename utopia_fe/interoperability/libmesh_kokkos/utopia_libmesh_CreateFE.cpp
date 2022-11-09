#include "utopia_libmesh_CreateFE.hpp"

#include "utopia_libmesh_RetroCompatibility.hpp"

// Dofs
#include <libmesh/dof_map.h>

// Mesh
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>

// FE
#include <libmesh/fe.h>
#include <libmesh/fe_base.h>
#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>

namespace utopia {

    void CreateFE<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        utopia::kokkos::FE<LibMeshScalar_t> &fe,
        const int order,
        const int var) {
        //
        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        if (ebegin == eend) {
            fe.clear();
            return;
        }

        auto elem_ptr = *ebegin;
        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(elem_ptr, dof_indices, var);

        int n_fun_x_elem = dof_indices.size();
        int spatial_dim = space.mesh().spatial_dimension();
        int mesh_dim = space.mesh().manifold_dimension();

        auto fe_type = dof_map.variable(var).type();
        auto quad = std::make_shared<libMesh::QGauss>(mesh_dim, libMesh::Order(order));

        std::unique_ptr<::libMesh::FEBase> fe_values = ::libMesh::FEBase::build(mesh_dim, fe_type);
        fe_values->attach_quadrature_rule(quad.get());

        auto &fun = fe_values->get_phi();
        auto &measure = fe_values->get_JxW();
        auto &grad = fe_values->get_dphi();

        fe_values->reinit(elem_ptr);

        int n_qp = measure.size();

        // Device buffers
        FE::MeasureView d_measure = FE::MeasureView("FE::MeasureView", n_cells, n_qp);
        FE::FunctionView d_fun = FE::FunctionView("FE::FunctionView", n_fun_x_elem, n_qp);
        FE::GradientView d_grad = FE::GradientView("FE::GradientView", n_cells, n_fun_x_elem, n_qp, spatial_dim);
        FE::IntView d_element_tags = FE::IntView("FE::IntView", n_cells);

        // Host mirror
        FE::MeasureView::HostMirror h_measure = Kokkos::create_mirror_view(d_measure);
        FE::FunctionView::HostMirror h_fun = Kokkos::create_mirror_view(d_fun);
        FE::GradientView::HostMirror h_grad = Kokkos::create_mirror_view(d_grad);
        FE::IntView::HostMirror h_element_tags = Kokkos::create_mirror_view(d_element_tags);

        // Shape functions are the same on all elements!
        for (int i = 0; i < n_fun_x_elem; ++i) {
            for (int q = 0; q < n_qp; ++q) {
                h_fun(i, q) = fun[i][q];
            }
        }

        SizeType eidx = -1;
        for (auto it = ebegin; it != eend; ++it) {
            auto elem_ptr = *it;
            ++eidx;

            int block = elem_ptr->subdomain_id();
            fe_values->reinit(elem_ptr);

            assert(n_qp == int(measure.size()));
            assert(n_qp == int(fun[0].size()));
            assert(n_qp == int(grad[0].size()));

            assert(n_fun_x_elem == int(fun.size()));
            assert(n_fun_x_elem == int(grad.size()));

            for (int q = 0; q < n_qp; ++q) {
                h_measure(eidx, q) = measure[q];
            }

            for (int i = 0; i < n_fun_x_elem; ++i) {
                for (int q = 0; q < n_qp; ++q) {
                    for (int d = 0; d < spatial_dim; ++d) {
                        h_grad(eidx, i, q, d) = grad[i][q](d);
                    }
                }
            }

            h_element_tags(eidx) = block;
        }

        // Move data to device
        Kokkos::deep_copy(d_measure, h_measure);
        Kokkos::deep_copy(d_fun, h_fun);
        Kokkos::deep_copy(d_grad, h_grad);
        Kokkos::deep_copy(d_element_tags, h_element_tags);

        fe.init(d_measure, d_fun, d_grad);
        fe.element_tags() = d_element_tags;
    }

    void get_surf_elem_dof_indices(const int sys_num,
                                   const int var_num,
                                   const int comp,
                                   const ::libMesh::Elem *elem_ptr,
                                   ::libMesh::Elem *side_ptr,
                                   std::vector<::libMesh::dof_id_type> &dof_indices) {
        const std::size_t n_side_nodes = side_ptr->n_nodes();
        const std::size_t n_vol_nodes = elem_ptr->n_nodes();

        dof_indices.clear();

        for (std::size_t k = 0; k < n_side_nodes; ++k) {
            const auto &surf_node = side_ptr->node_ref(k);
            if (!surf_node.has_dofs()) break;

            for (std::size_t j = 0; j < n_vol_nodes; ++j) {
                const auto &vol_node = elem_ptr->node_ref(j);

                if (vol_node.absolute_fuzzy_equals(surf_node, 1e-14)) {
                    const auto v_dof = vol_node.dof_number(sys_num, var_num, comp);

                    dof_indices.push_back(v_dof);
                    break;
                }
            }
        }
    }

    void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        FE &fe,
        int order,
        int var) {
        // All!
        apply(space, fe, "", order, var);
    }

    void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::kokkos::FE<LibMeshScalar_t>>::apply(
        const utopia::libmesh::FunctionSpace &space,
        FE &fe,
        const std::string &boundary_name,
        int order,
        int var) {
        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const int sys_num = space.system_id();
        const int comp = 0;

        ::libMesh::FEType fe_type = dof_map.variable_type(var);
        int n_cells = 0;
        int spatial_dim = space.mesh().spatial_dimension();
        int mesh_dim = space.mesh().manifold_dimension();

        auto &&binfo = mesh.get_boundary_info();

        int boundary_id_select = -1;
        if (!boundary_name.empty()) {
            boundary_id_select = binfo.get_id_by_name(boundary_name);
        }

        libMesh::Elem *boundary_element = nullptr;
        int boundary_side = 0;

        for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
            const int n_sides = elem_ptr->n_sides();

            for (int i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                bool selected =
                    boundary_id_select == -1 || boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr, i);

                if (!selected) continue;

                n_cells++;

                if (!boundary_element) {
                    // We have found a boundary element!
                    boundary_element = elem_ptr;
                    boundary_side = i;
                }
            }
        }

        if (n_cells == 0) {
            // Nothing to do here!
            fe.clear();
            return;
        }

        std::unique_ptr<libMesh::Elem> side_ptr = boundary_element->build_side_ptr(boundary_side);

        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(side_ptr.get(), dof_indices, var);

        int n_fun_x_elem = dof_indices.size();

        auto quad = std::make_shared<libMesh::QGauss>(mesh_dim - 1, libMesh::Order(order));

        std::unique_ptr<::libMesh::FEBase> fe_values = ::libMesh::FEBase::build(mesh_dim - 1, fe_type);
        fe_values->attach_quadrature_rule(quad.get());

        auto &fun = fe_values->get_phi();
        auto &measure = fe_values->get_JxW();
        auto &grad = fe_values->get_dphi();

        fe_values->reinit(side_ptr.get());

        int n_qp = measure.size();

        // // Device buffers
        FE::MeasureView d_measure = FE::MeasureView("FE::MeasureView", n_cells, n_qp);
        FE::FunctionView d_fun = FE::FunctionView("FE::FunctionView", n_fun_x_elem, n_qp);
        FE::GradientView d_grad = FE::GradientView("FE::GradientView", n_cells, n_fun_x_elem, n_qp, spatial_dim);
        FE::IntView d_element_tags = FE::IntView("FE::IntView", n_cells);

        // // Host mirror
        FE::MeasureView::HostMirror h_measure = Kokkos::create_mirror_view(d_measure);
        FE::FunctionView::HostMirror h_fun = Kokkos::create_mirror_view(d_fun);
        FE::GradientView::HostMirror h_grad = Kokkos::create_mirror_view(d_grad);
        FE::IntView::HostMirror h_element_tags = Kokkos::create_mirror_view(d_element_tags);

        // double totm = 0.0;

        int eidx = 0;
        for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
            const int n_sides = elem_ptr->n_sides();

            for (int i = 0; i < n_sides; ++i) {
                if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                    continue;
                }

                bool selected =
                    boundary_id_select == -1 || boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr, i);

                if (!selected) continue;

                int sideset = utopia::libmesh::boundary_id(binfo, elem_ptr, i);

                std::unique_ptr<libMesh::Elem> side_ptr = elem_ptr->build_side_ptr(i);
                fe_values->reinit(side_ptr.get());

                assert(n_qp == int(measure.size()));
                assert(n_qp == int(fun[0].size()));
                assert(n_qp == int(grad[0].size()));

                assert(n_fun_x_elem == int(fun.size()));
                assert(n_fun_x_elem == int(grad.size()));

                for (int q = 0; q < n_qp; ++q) {
                    h_measure(eidx, q) = measure[q];
                }

                for (int i = 0; i < n_fun_x_elem; ++i) {
                    for (int q = 0; q < n_qp; ++q) {
                        for (int d = 0; d < spatial_dim; ++d) {
                            h_grad(eidx, i, q, d) = grad[i][q](d);
                        }
                    }
                }

                h_element_tags(eidx) = sideset;
                eidx++;
            }
        }

        // Move data to device
        Kokkos::deep_copy(d_measure, h_measure);
        Kokkos::deep_copy(d_fun, h_fun);
        Kokkos::deep_copy(d_grad, h_grad);
        Kokkos::deep_copy(d_element_tags, h_element_tags);

        fe.init(d_measure, d_fun, d_grad);
        fe.element_tags() = d_element_tags;
    }

    void ConvertField<utopia::Field<utopia::libmesh::FunctionSpace>, LibMeshFEField_t>::apply(const FromField &from,
                                                                                              ToField &to) {
        assert(false);
        Utopia::Abort();
    }

    void ConvertField<LibMeshFEField_t, utopia::Field<utopia::libmesh::FunctionSpace>>::apply(const FromField &from,
                                                                                              ToField &to) {
        assert(false);
        Utopia::Abort();
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscMatrix>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix,
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscMatrix>::Scalar;
        using IndexArray = Traits<PetscMatrix>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        // Scoped lock for parallel assembly
        Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*ebegin, dof_indices, var);

        IndexArray idx(dof_indices.size(), 0);
        std::vector<Scalar> values(dof_indices.size() * dof_indices.size(), 0);

        LibMeshViewDevice_t::HostMirror h_element_matrices = Kokkos::create_mirror_view(d_element_matrices);
        Kokkos::deep_copy(h_element_matrices, d_element_matrices);

        int nfun_x_nvar = h_element_matrices.extent(1);
        int n_var = nfun_x_nvar / dof_indices.size();

        assert(var + n_var <= space.n_var());

        SizeType eidx = -1;
        for (auto it = ebegin; it != eend; ++it) {
            auto elem_ptr = *it;
            ++eidx;

            for (int v = 0; v < n_var; ++v) {
                dof_map.dof_indices(*ebegin, dof_indices, var + v);

                int n_shape_functions = dof_indices.size();

                // Convert to PETSc types
                for (int i = 0; i < n_shape_functions; ++i) {
                    idx[i] = dof_indices[i];
                }

                for (int i = 0; i < n_shape_functions; ++i) {
                    for (int j = 0; j < n_shape_functions; ++j) {
                        auto value = h_element_matrices(eidx, i * n_var + v, j * n_var + v);
                        values[i * n_shape_functions + j] = value;
                    }
                }

                matrix.add_matrix(idx, idx, values);
            }
        }
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        int var) {
        apply(space, d_element_vectors, mode, vector, -1, var);
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int check_n_var,  // n_var is automatically determined inside!
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscVector>::Scalar;
        using IndexArray = Traits<PetscVector>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        // Scoped lock for parallel assembly
        Write<PetscVector> w(vector, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        std::vector<::libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*ebegin, dof_indices, var);

        IndexArray idx(dof_indices.size(), 0);
        std::vector<Scalar> values(dof_indices.size(), 0);

        LibMeshViewDevice_t::HostMirror h_element_vectors = Kokkos::create_mirror_view(d_element_vectors);
        Kokkos::deep_copy(h_element_vectors, d_element_vectors);

        int nfun_x_nvar = h_element_vectors.extent(1);
        int n_var = nfun_x_nvar / dof_indices.size();

        assert(check_n_var == -1 || n_var == check_n_var);

        assert(var + n_var <= space.n_var());

        SizeType eidx = -1;
        for (auto it = ebegin; it != eend; ++it) {
            auto elem_ptr = *it;
            ++eidx;

            for (int v = 0; v < n_var; ++v) {
                dof_map.dof_indices(elem_ptr, dof_indices, var + v);

                int n_shape_functions = dof_indices.size();

                // Convert to PETSc types
                for (int i = 0; i < n_shape_functions; ++i) {
                    idx[i] = dof_indices[i];
                }

                for (int i = 0; i < n_shape_functions; ++i) {
                    auto value = h_element_vectors(eidx, i * n_var + v);
                    values[i] = value;
                }

                vector.add_vector(idx, values);
            }
        }
    }

    void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t, PetscVector>::side_apply(
        const utopia::libmesh::FunctionSpace &space,
        const LibMeshViewDevice_t &d_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &boundary_name,
        int var) {
        // Petsc type defs
        using Scalar = Traits<PetscVector>::Scalar;
        using IndexArray = Traits<PetscVector>::IndexArray;

        auto &&mesh = space.mesh().raw_type();
        auto &&dof_map = space.raw_type_dof_map();

        const auto n_local_dofs = dof_map.n_local_dofs();
        const long n_cells = mesh.n_active_local_elem();
        auto ebegin = mesh.active_local_elements_begin();
        auto eend = mesh.active_local_elements_end();

        auto &&binfo = mesh.get_boundary_info();

        int boundary_id_select = -1;
        if (!boundary_name.empty()) {
            boundary_id_select = binfo.get_id_by_name(boundary_name);
        }

        // Scoped lock for parallel assembly
        Write<PetscVector> w(vector, utopia::GLOBAL_ADD);
        if (ebegin == eend) return;

        libMesh::Elem *boundary_element = nullptr;
        int boundary_side = 0;

        {
            // Find one boundary element!
            for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
                const int n_sides = elem_ptr->n_sides();

                for (int i = 0; i < n_sides; ++i) {
                    if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                        continue;
                    }

                    bool selected = boundary_id_select == -1 ||
                                    boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr, i);

                    if (!selected) continue;

                    if (!boundary_element) {
                        // We have found a boundary element!
                        boundary_element = elem_ptr;
                        boundary_side = i;
                        break;
                    }
                }

                if (boundary_element) break;
            }
        }

        if (!boundary_element) return;

        // std::unique_ptr<libMesh::Elem> side_ptr = boundary_element->build_side_ptr(boundary_side);

        // std::vector<::libMesh::dof_id_type> dof_indices;
        // dof_map.dof_indices(side_ptr, dof_indices, var);

        // IndexArray idx(dof_indices.size(), 0);
        // std::vector<Scalar> values(dof_indices.size(), 0);

        // LibMeshViewDevice_t::HostMirror h_element_vectors = Kokkos::create_mirror_view(d_element_vectors);
        // Kokkos::deep_copy(h_element_vectors, d_element_vectors);

        // int nfun_x_nvar = h_element_vectors.extent(1);
        // int n_var = nfun_x_nvar / dof_indices.size();

        // assert(var + n_var <= space.n_var());

        // SizeType eidx = 0;
        // for (auto it = ebegin; it != eend; ++it) {
        //     auto elem_ptr = *it;

        //     for (int i = 0; i < n_sides; ++i) {
        //         if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
        //             continue;
        //         }

        //         bool selected =
        //             boundary_id_select == -1 || boundary_id_select == utopia::libmesh::boundary_id(binfo, elem_ptr,
        //             i);

        //         if (!selected) continue;

        //         std::unique_ptr<libMesh::Elem> side_ptr = elem_ptr->build_side_ptr(i);

        //         for (int v = 0; v < n_var; ++v) {
        //             dof_map.dof_indices(side_ptr.get(), dof_indices, var + v);

        //             int n_shape_functions = dof_indices.size();

        //             // Convert to PETSc types
        //             for (int i = 0; i < n_shape_functions; ++i) {
        //                 idx[i] = dof_indices[i];
        //             }

        //             for (int i = 0; i < n_shape_functions; ++i) {
        //                 auto value = h_element_vectors(eidx, i * n_var + v);
        //                 values[i] = value;
        //             }

        //             vector.add_vector(idx, values);
        //             eidx++;
        //         }
        //     }
        // }
    }

}  // namespace utopia
