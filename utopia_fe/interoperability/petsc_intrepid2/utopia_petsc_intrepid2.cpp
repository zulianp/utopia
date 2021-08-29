#include "utopia_petsc_intrepid2.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>

#include "utopia_intrepid2_FE.hpp"

namespace utopia {

    static ::shards::CellTopology convert_elem_type(int dim, int n_elements_cell) {
        switch (dim) {
            case 2: {
                switch (n_elements_cell) {
                    case 1: {
                        return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
                    }
                    case 2: {
                        return ::shards::getCellTopologyData<::shards::Triangle<>>();
                    }
                    default: {
                        break;
                    }
                }
            }
            case 3: {
                switch (n_elements_cell) {
                    case 1: {
                        return ::shards::getCellTopologyData<::shards::Hexahedron<>>();
                    }
                    case 2: {
                        return ::shards::getCellTopologyData<::shards::Tetrahedron<>>();
                    }
                    default: {
                        break;
                    }
                }
            }
            default: {
                assert(false);
                return ::shards::getCellTopologyData<::shards::Node>();
            }
        }
    }

    template <typename Scalar>
    void CreateFE<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::petsc::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        auto &mesh = space.mesh();
        auto mesh_view = mesh.view();
        SizeType n_local_elements = mesh_view.n_local_elements();
        int spatial_dim = mesh_view.spatial_dim();
        int n_nodes_x_elem = mesh_view.n_nodes_x_element();
        auto topo = convert_elem_type(mesh_view.dim(), mesh_view.elements_x_cell());

        PetscViewDevice_t<Scalar> device_cell_points("cell_points", n_local_elements, n_nodes_x_elem, spatial_dim);
        PetscIntViewDevice_t device_element_tags("element_tags", n_local_elements);

        typename PetscViewDevice_t<Scalar>::HostMirror cell_points = Kokkos::create_mirror_view(device_cell_points);
        PetscIntViewDevice_t::HostMirror element_tags = Kokkos::create_mirror_view(device_element_tags);

        auto elems = mesh.elem_to_node_index();
        // auto local_elems = mesh.elem_to_local_node_index();
        // TODO Create points and tags

        utopia::petsc::StructuredGrid::View::Point p(spatial_dim);

        for (SizeType idx = 0; idx < n_local_elements; ++idx) {
            for (int k = 0; k < n_nodes_x_elem; ++k) {
                auto node_idx = elems(idx, k);
                mesh_view.point(node_idx, p);

                // utopia::out() << node_idx << ") ";

                // for (auto pi : p) {
                //     utopia::out() << pi << " ";
                // }

                // utopia::out() << '\n';

                for (int d = 0; d < spatial_dim; ++d) {
                    cell_points(idx, k, d) = p[d];
                }
            }
        }

        Kokkos::deep_copy(device_cell_points, cell_points);
        Kokkos::deep_copy(device_element_tags, element_tags);

        fe.init(topo, device_cell_points, degree);
        fe.element_tags() = device_element_tags;
    }

    // Explicit instantation
    template class CreateFE<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<PetscScalar_t>>;

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::petsc::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        assert(false);
        Utopia::Abort("IMPLEMENT ME!");
    }

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::petsc::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const std::string &part_name,
        const int degree) {
        assert(false);
        Utopia::Abort("IMPLEMENT ME!");
    }

    template class CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<PetscScalar_t>>;

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscMatrix>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &device_element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix) {
        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2)");

        // // Type defs
        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
        using Size_t = Traits<PetscMatrix>::SizeType;

        auto &mesh = space.mesh();
        auto mesh_view = mesh.view();

        SizeType n_local_elements = mesh_view.n_local_elements();
        int spatial_dim = mesh_view.spatial_dim();
        int n_nodes_x_elem = mesh_view.n_nodes_x_element();

        typename PetscViewDevice_t<Scalar>::HostMirror element_matrices =
            ::Kokkos::create_mirror_view(device_element_matrices);
        ::Kokkos::deep_copy(element_matrices, device_element_matrices);

        const int n_var = space.n_var();

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        if (matrix.empty()) {
            space.create_matrix(matrix);
        } else {
            // Reuse matrix
            if (matrix.is_assembled() && mode == OVERWRITE_MODE) {
                matrix *= 0.0;
            } else if (mode == SUBTRACT_MODE) {
                matrix *= -1.0;
            }
        }

        bool is_block = matrix.is_block();

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");
        {
            Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

            const SizeType n_dofs = element_matrices.extent(1);
            const SizeType nn = n_dofs / n_var;

            IndexArray_t idx(is_block ? nn : n_dofs);
            ScalarArray_t val(n_dofs * n_dofs);

            const bool is_block = matrix.is_block();

            auto elems = mesh.elem_to_node_index();

            for (SizeType elem_idx = 0; elem_idx < n_local_elements; ++elem_idx) {
                for (Size_t i = 0, dof_k = 0; i < nn; ++i) {
                    if (!is_block) {
                        for (int v = 0; v < n_var; ++v, ++dof_k) {
                            idx[dof_k] = elems(elem_idx, i) * n_var + v;
                        }
                    } else {
                        idx[i] = elems(elem_idx, i);
                    }
                }

                for (Size_t i = 0; i < nn; ++i) {
                    for (int di = 0; di < n_var; ++di) {
                        const Size_t idx_i = i * n_var + di;

                        for (Size_t j = 0; j < nn; ++j) {
                            for (int dj = 0; dj < n_var; ++dj) {
                                const Size_t idx_j = j * n_var + dj;

                                val[idx_i * n_dofs + idx_j] = element_matrices(elem_idx, idx_i, idx_j);
                            }
                        }
                    }
                }

                if (is_block) {
                    matrix.add_matrix_blocked(idx, idx, val);
                } else {
                    matrix.add_matrix(idx, idx, val);
                }
            }
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");

        if (mode == SUBTRACT_MODE) {
            matrix *= -1.0;
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2)");
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector) {
        LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::apply(
            space, element_vectors, mode, vector, space.n_var());
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &device_element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int n_var) {
        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2)");

        // // Type defs
        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
        using Size_t = Traits<PetscMatrix>::SizeType;

        auto &mesh = space.mesh();
        auto mesh_view = mesh.view();

        SizeType n_local_elements = mesh_view.n_local_elements();
        int spatial_dim = mesh_view.spatial_dim();
        int n_nodes_x_elem = mesh_view.n_nodes_x_element();

        typename PetscViewDevice_t<Scalar>::HostMirror element_vectors =
            ::Kokkos::create_mirror_view(device_element_vectors);
        ::Kokkos::deep_copy(element_vectors, device_element_vectors);

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        if (vector.empty()) {
            space.create_vector(vector);
        } else {
            // Reuse vector
            if (mode == OVERWRITE_MODE) {
                vector *= 0.0;
            } else if (mode == SUBTRACT_MODE) {
                vector *= -1.0;
            }
        }

        vector.set_block_size(n_var);

        const bool is_block = vector.block_size() > 1;

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");
        {
            Write<PetscVector> w(vector, utopia::GLOBAL_ADD);

            const SizeType n_dofs = element_vectors.extent(1);
            const SizeType nn = n_dofs / n_var;

            IndexArray_t idx(is_block ? nn : n_dofs);
            ScalarArray_t val(n_dofs * n_dofs);

            auto elems = mesh.elem_to_node_index();

            for (SizeType elem_idx = 0; elem_idx < n_local_elements; ++elem_idx) {
                for (Size_t i = 0, dof_k = 0; i < nn; ++i) {
                    if (!is_block) {
                        for (int v = 0; v < n_var; ++v, ++dof_k) {
                            idx[dof_k] = elems(elem_idx, i) * n_var + v;
                        }
                    } else {
                        idx[i] = elems(elem_idx, i);
                    }
                }

                for (Size_t i = 0; i < nn; ++i) {
                    for (int di = 0; di < n_var; ++di) {
                        const Size_t idx_i = i * n_var + di;
                        val[idx_i] = element_vectors(elem_idx, idx_i);
                    }
                }

                if (is_block) {
                    vector.add_vector_blocked(idx, val);
                } else {
                    vector.add_vector(idx, val);
                }
            }
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");

        if (mode == SUBTRACT_MODE) {
            vector *= -1.0;
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2)");
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::side_apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &part_name) {
        assert(false);
        Utopia::Abort("IMPLEMENT ME!");
    }

    template class LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<PetscScalar_t>, PetscMatrix>;
    template class LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<PetscScalar_t>, PetscVector>;

    template <typename Scalar>
    void ConvertField<Field<utopia::petsc::FunctionSpace>, Intrepid2Field<Scalar>>::apply(
        const Field<utopia::petsc::FunctionSpace> &from,
        Intrepid2Field<Scalar> &to) {
        GlobalToLocal<utopia::petsc::FunctionSpace, Vector, PetscViewDevice_t<Scalar>>::apply(
            *from.space(), from.data(), to.data(), from.tensor_size());

        to.set_tensor_size(from.tensor_size());
    }

    template class ConvertField<Field<utopia::petsc::FunctionSpace>, Intrepid2Field<PetscScalar_t>>;

    template <typename Scalar>
    void GlobalToLocal<utopia::petsc::FunctionSpace,
                       Traits<utopia::petsc::FunctionSpace>::Vector,
                       PetscViewDevice_t<Scalar>>::apply(const utopia::petsc::FunctionSpace &space,
                                                         const Vector &vector,
                                                         PetscViewDevice_t<Scalar> &device_element_vectors,
                                                         const int n_comp) {
        UTOPIA_TRACE_REGION_BEGIN("GlobalToLocal(Petsc,Intrepid2)");

        Vector local;
        space.global_to_local(vector, local);

        auto &mesh = space.mesh();
        auto mesh_view = mesh.view();
        SizeType n_local_elements = mesh_view.n_local_elements();
        int spatial_dim = mesh_view.spatial_dim();
        int n_nodes_x_elem = mesh_view.n_nodes_x_element();
        int n_var = space.n_var();
        auto topo = convert_elem_type(mesh_view.dim(), mesh_view.elements_x_cell());

        if (Size_t(device_element_vectors.extent(0)) < n_local_elements ||
            Size_t(device_element_vectors.extent(1)) < n_nodes_x_elem * n_comp) {
            device_element_vectors =
                PetscViewDevice_t<Scalar>("Coefficients", n_local_elements, n_nodes_x_elem * n_comp);
        }

        typename PetscViewDevice_t<Scalar>::HostMirror element_vectors =
            ::Kokkos::create_mirror_view(device_element_vectors);

        // FIXME not on device
        auto local_view = const_local_view_device(local);

        auto local_elems = mesh.elem_to_local_node_index();

        for (SizeType idx = 0; idx < n_local_elements; ++idx) {
            for (int k = 0; k < n_nodes_x_elem; ++k) {
                auto node_idx = local_elems(idx, k);

                for (int i = 0; i < n_nodes_x_elem; ++i) {
                    for (int d = 0; d < n_comp; ++d) {
                        element_vectors(idx, i * n_comp + d) = local_view.get(node_idx * n_comp + d);
                    }
                }
            }
        }

        ::Kokkos::deep_copy(device_element_vectors, element_vectors);

        UTOPIA_TRACE_REGION_END("GlobalToLocal(Petsc,Intrepid2)");
    }

    template class GlobalToLocal<utopia::petsc::FunctionSpace,
                                 Traits<utopia::petsc::FunctionSpace>::Vector,
                                 PetscViewDevice_t<PetscScalar_t>>;

    void l2_norm(const Field<utopia::petsc::FunctionSpace> &field, std::vector<PetscScalar_t> &norms) {
        using Scalar = typename Traits<utopia::petsc::FunctionSpace>::Scalar;

        assert(false);
        Utopia::Abort("IMPLEMENT ME!");
        // utopia::kokkos::L2Norm<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>().compute(field, norms);
    }

}  // namespace utopia