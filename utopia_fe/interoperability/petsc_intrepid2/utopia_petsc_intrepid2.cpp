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
        auto local_elems = mesh.elem_to_local_node_index();
        // TODO Create points and tags

        utopia::petsc::StructuredGrid::View::Point p(spatial_dim);

        for (SizeType idx = 0; idx < n_local_elements; ++idx) {
            for (int k = 0; k < n_nodes_x_elem; ++k) {
                auto node_idx = elems(idx, k);
                mesh_view.point(idx, p);

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
        const int degree) {}

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::petsc::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const std::string &part_name,
        const int degree) {}

    template class CreateFEOnBoundary<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<PetscScalar_t>>;

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscMatrix>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &device_element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix) {
        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2)");

        // // Type defs
        // using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        // using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
        // using Size_t = Traits<PetscMatrix>::SizeType;

        // using Bucket_t = ::petsc::mesh::Bucket;
        // using BucketVector_t = ::petsc::mesh::BucketVector;
        // using Entity_t = ::petsc::mesh::Entity;

        // typename PetscViewDevice_t<Scalar>::HostMirror element_matrices =
        //     ::Kokkos::create_mirror_view(device_element_matrices);
        // ::Kokkos::deep_copy(element_matrices, device_element_matrices);

        // const int n_var = space.n_var();

        // UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        // if (!matrix.is_block() && n_var != 1) {
        //     matrix.clear();
        // }

        // if (matrix.empty()) {
        //     space.create_matrix(matrix);

        // } else {
        //     // Reuse matrix
        //     if (matrix.is_assembled() && mode == OVERWRITE_MODE) {
        //         matrix *= 0.0;
        //     } else if (mode == SUBTRACT_MODE) {
        //         matrix *= -1.0;
        //     }
        // }

        // // Fix preallocation bug and REMOVE ME
        // if (!space.mesh().has_aura() && space.n_var() > 1 && space.comm().size() > 1) {
        //     MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        //     m_utopia_warning(
        //         "using: MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR,"
        //         "PETSC_FALSE);");
        // }

        // UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix handling)");

        // auto &bulk_data = space.mesh().bulk_data();

        // const BucketVector_t &elem_buckets = utopia::petsc::local_elements(bulk_data);

        // UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");
        // {
        //     Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

        //     const SizeType n_dofs = element_matrices.extent(1);
        //     const SizeType nn = n_dofs / n_var;

        //     IndexArray_t idx(nn);
        //     ScalarArray_t val(n_dofs * n_dofs);

        //     const bool is_block = matrix.is_block();

        //     auto &&local_to_global = space.dof_map().local_to_global();

        //     Size_t elem_idx = 0;
        //     for (const auto &ib : elem_buckets) {
        //         const Bucket_t &b = *ib;
        //         const Bucket_t::size_type length = b.size();

        //         for (Bucket_t::size_type k = 0; k < length; ++k) {
        //             // get the current node entity and extract the id to fill it into the field
        //             Entity_t elem = b[k];
        //             // const Size_t elem_idx = utopia::petsc::convert_entity_to_index(elem) - n_local_nodes;
        //             const Size_t n_nodes = bulk_data.num_nodes(elem);
        //             UTOPIA_UNUSED(n_nodes);

        //             assert(nn == n_nodes);

        //             auto node_ids = bulk_data.begin_nodes(elem);

        //             if (local_to_global.empty()) {
        //                 for (Size_t i = 0; i < nn; ++i) {
        //                     idx[i] = utopia::petsc::convert_entity_to_index(node_ids[i]);
        //                     assert(idx[i] < space.n_dofs());
        //                     assert(idx[i] >= 0);
        //                 }
        //             } else {
        //                 for (Size_t i = 0; i < nn; ++i) {
        //                     idx[i] = local_to_global.block(utopia::petsc::convert_entity_to_index(node_ids[i]));
        //                     assert(idx[i] < space.n_dofs());
        //                     assert(idx[i] >= 0);
        //                 }
        //             }

        //             for (Size_t i = 0; i < nn; ++i) {
        //                 for (int di = 0; di < n_var; ++di) {
        //                     const Size_t idx_i = i * n_var + di;

        //                     for (Size_t j = 0; j < nn; ++j) {
        //                         for (int dj = 0; dj < n_var; ++dj) {
        //                             const Size_t idx_j = j * n_var + dj;

        //                             val[idx_i * n_dofs + idx_j] = element_matrices(elem_idx, idx_i, idx_j);
        //                         }
        //                     }
        //                 }
        //             }

        //             if (is_block) {
        //                 matrix.add_matrix_blocked(idx, idx, val);
        //             } else {
        //                 matrix.add_matrix(idx, idx, val);
        //             }

        //             ++elem_idx;
        //         }
        //     }
        // }
        // UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2,Matrix assembly)");

        // if (mode == SUBTRACT_MODE) {
        //     matrix *= -1.0;
        // }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Petsc,Intrepid2)");
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector) {
        // LocalToGlobalFromBuckets<Scalar>::apply(
        // space, utopia::petsc::local_elements(space.mesh().bulk_data()), element_vectors, mode, vector);
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int n_var) {
        // LocalToGlobalFromBuckets<Scalar>::apply(
        // space, utopia::petsc::local_elements(space.mesh().bulk_data()), element_vectors, mode, vector, n_var);
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::petsc::FunctionSpace, PetscViewDevice_t<Scalar>, PetscVector>::side_apply(
        const utopia::petsc::FunctionSpace &space,
        const PetscViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &part_name) {}

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

        // auto &bulk_data = space.mesh().bulk_data();
        // const auto &elem_buckets = utopia::petsc::local_elements(bulk_data);
        // const Size_t num_elem = utopia::petsc::count_local_elements(bulk_data);

        // // Dirty hack (FIXME once petsc usage is a bit more profficient)
        // auto *first_bucket = *elem_buckets.begin();
        // auto topo = convert_elem_type(first_bucket->topology(), true);
        // Size_t n_nodes_x_elem = bulk_data.num_nodes((*first_bucket)[0]);

        // if (Size_t(device_element_vectors.extent(0)) < num_elem ||
        //     Size_t(device_element_vectors.extent(1)) < n_nodes_x_elem * n_comp) {
        //     device_element_vectors = PetscViewDevice_t<Scalar>("Coefficients", num_elem, n_nodes_x_elem * n_comp);
        // }

        typename PetscViewDevice_t<Scalar>::HostMirror element_vectors =
            ::Kokkos::create_mirror_view(device_element_vectors);

        // FIXME not on device
        auto local_view = const_local_view_device(local);

        // Size_t elem_idx = 0;
        // for (const auto &ib : elem_buckets) {
        //     const auto &b = *ib;
        //     const SizeType length = b.size();

        //     for (SizeType k = 0; k < length; ++k) {
        //         // get the current node entity and extract the id to fill it into the field
        //         auto elem = b[k];
        //         const int n_nodes = bulk_data.num_nodes(elem);

        //         auto node_ids = bulk_data.begin_nodes(elem);

        //         for (int i = 0; i < n_nodes; ++i) {
        //             SizeType node_idx = utopia::petsc::convert_entity_to_index(node_ids[i]);
        //             for (int d = 0; d < n_comp; ++d) {
        //                 element_vectors(elem_idx, i * n_comp + d) = local_view.get(node_idx * n_comp + d);
        //             }
        //         }

        //         ++elem_idx;
        //     }
        // }

        ::Kokkos::deep_copy(device_element_vectors, element_vectors);

        UTOPIA_TRACE_REGION_END("GlobalToLocal(Petsc,Intrepid2)");
    }

    template class GlobalToLocal<utopia::petsc::FunctionSpace,
                                 Traits<utopia::petsc::FunctionSpace>::Vector,
                                 PetscViewDevice_t<PetscScalar_t>>;

    void l2_norm(const Field<utopia::petsc::FunctionSpace> &field, std::vector<PetscScalar_t> &norms) {
        using Scalar = typename Traits<utopia::petsc::FunctionSpace>::Scalar;
        // utopia::kokkos::L2Norm<utopia::petsc::FunctionSpace, utopia::intrepid2::FE<Scalar>>().compute(field, norms);
    }

}  // namespace utopia