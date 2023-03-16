// #include "utopia_libmesh_intrepid2.hpp"

// #include "utopia_libmesh_intrepid2.hpp"
// #include "utopia_petsc_Matrix.hpp"

// #include "utopia_intrepid2_FE.hpp"
// #include "utopia_libmesh_FunctionSpace_new.hpp"

// #include "libmesh/dof_map.h"
// #include "libmesh/elem.h"
// #include "libmesh/mesh_base.h"

// namespace utopia {

//     static ::shards::CellTopology convert_elem_type(const libMesh::ElemType &type) {
//         using ::shards::getCellTopologyData;

//         switch (type) {
//             case libMesh::EDGE2:
//                 return getCellTopologyData<::shards::Line<>>();
//             case libMesh::EDGE3:
//                 return getCellTopologyData<::shards::Line<3>>();
//             case libMesh::TRI3:
//                 return getCellTopologyData<::shards::Triangle<>>();
//             case libMesh::TRI6:
//                 return getCellTopologyData<::shards::Triangle<6>>();
//             case libMesh::QUAD4:
//                 return getCellTopologyData<::shards::Quadrilateral<>>();
//             case libMesh::QUAD8:
//                 return getCellTopologyData<::shards::Quadrilateral<8>>();
//             case libMesh::QUAD9:
//                 return getCellTopologyData<::shards::Quadrilateral<9>>();
//             case libMesh::TET4:
//                 return getCellTopologyData<::shards::Tetrahedron<>>();
//             case libMesh::TET10:
//                 return getCellTopologyData<::shards::Tetrahedron<10>>();
//             case libMesh::HEX8:
//                 return getCellTopologyData<::shards::Hexahedron<>>();
//             case libMesh::HEX27:
//                 return getCellTopologyData<::shards::Hexahedron<27>>();
//             case libMesh::QUADSHELL8:
//                 return getCellTopologyData<::shards::Quadrilateral<8>>();
//             case libMesh::QUADSHELL4:
//                 return getCellTopologyData<::shards::Quadrilateral<4>>();
//             default: {
//                 assert(false);
//                 return getCellTopologyData<::shards::Node>();
//             }
//         }
//     }

//     template <typename Scalar>
//     void CreateFE<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
//         const utopia::libmesh::FunctionSpace &space,
//         utopia::intrepid2::FE<Scalar> &fe,
//         const int degree) {
//         // using FE = utopia::intrepid2::FE<Scalar>;

//         auto &m = space.mesh().raw_type();

//         auto e_begin = m.active_local_elements_begin();
//         auto e_end = m.active_local_elements_end();
//         if (e_begin != e_end) {
//             auto n_elems = m.n_active_local_elem();
//             auto n_nodes_x_elem = (*m.active_local_elements_begin())->n_nodes();
//             auto topo = convert_elem_type((*m.active_local_elements_begin())->type());
//             int spatial_dim = m.spatial_dimension();

//             LibMeshViewDevice_t<Scalar> device_cell_points("cell_points", n_elems, n_nodes_x_elem, spatial_dim);
//             typename LibMeshViewDevice_t<Scalar>::HostMirror cell_points =
//                 ::Kokkos::create_mirror_view(device_cell_points);

//             LibMeshIntViewDevice_t device_element_tags("element_tags", n_elems);
//             LibMeshIntViewDevice_t::HostMirror element_tags = ::Kokkos::create_mirror_view(device_element_tags);

//             int index = 0;

//             for (auto e_it = e_begin; e_it != e_end; ++e_it, ++index) {
//                 const auto &e = **e_it;

//                 element_tags(index) = e.subdomain_id();
//                 const int n_nodes = e.n_nodes();
//                 for (int i = 0; i < n_nodes; ++i) {
//                     const auto &node = e.node_ref(i);
//                     for (int d = 0; d < spatial_dim; ++d) {
//                         cell_points(index, i, d) = node(d);
//                     }
//                 }
//             }

//             Kokkos::deep_copy(device_cell_points, cell_points);
//             Kokkos::deep_copy(device_element_tags, element_tags);

//             fe.init(topo, device_cell_points, degree);
//             fe.element_tags = device_element_tags;
//         }
//     }

//     // // Explicit instantation
//     template class CreateFE<utopia::libmesh::FunctionSpace,
//                             utopia::intrepid2::FE<Traits<utopia::libmesh::FunctionSpace>::Scalar>>;

//     template <typename Scalar>
//     void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
//         const utopia::libmesh::FunctionSpace &space,
//         utopia::intrepid2::FE<Scalar> &fe,
//         const int degree) {
//         assert(false);
//     }

//     template <typename Scalar>
//     void CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
//         const utopia::libmesh::FunctionSpace &space,
//         utopia::intrepid2::FE<Scalar> &fe,
//         const std::string &part_name,
//         const int degree) {
//         assert(false);
//     }

//     template class CreateFEOnBoundary<utopia::libmesh::FunctionSpace, utopia::intrepid2::FE<double>>;

//     template <typename Scalar>
//     void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<Scalar>, PetscMatrix>::apply(
//         const utopia::libmesh::FunctionSpace &space,
//         const LibMeshViewDevice_t<Scalar> &device_element_matrices,
//         AssemblyMode mode,
//         PetscMatrix &matrix) {
//         UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(libMesh,Intrepid2)");

//         // Type defs
//         using IndexArray_t = Traits<PetscMatrix>::IndexArray;
//         using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
//         using Size_t = Traits<PetscMatrix>::SizeType;
//         using HostView = typename LibMeshViewDevice_t<Scalar>::HostMirror;

//         HostView element_matrices = ::Kokkos::create_mirror_view(device_element_matrices);
//         ::Kokkos::deep_copy(element_matrices, device_element_matrices);

//         if (matrix.empty()) {
//             space.create_matrix(matrix);

//         } else {
//             // Reuse matrix
//             if (matrix.is_assembled() && mode == OVERWRITE_MODE) {
//                 matrix *= 0.0;
//             } else if (mode == SUBTRACT_MODE) {
//                 matrix *= -1.0;
//             }
//         }

//         // Fix preallocation bug and REMOVE ME
//         // if (space.n_var() > 1 && space.comm().size() > 1) {
//         //     MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//         //     m_utopia_warning(
//         //         "using: MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR,"
//         //         "PETSC_FALSE);");
//         // }

//         auto &m = space.mesh().raw_type();
//         auto &dof_map = space.raw_type_dof_map();

//         {
//             Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

//             auto e_begin = m.active_local_elements_begin();
//             auto e_end = m.active_local_elements_end();

//             if (e_begin != e_end) {
//                 const int n_var = space.n_var();
//                 const SizeType n_dofs = element_matrices.extent(1);
//                 const SizeType nn = n_dofs / n_var;

//                 const bool is_block = matrix.is_block();
//                 IndexArray_t idx(is_block ? nn : n_dofs);
//                 ScalarArray_t val(n_dofs * n_dofs);
//                 std::vector<libMesh::dof_id_type> indices;

//                 int elem_idx = 0;
//                 for (auto e_it = e_begin; e_it != e_end; ++e_it, ++elem_idx) {
//                     const auto &e = **e_it;

//                     dof_map.dof_indices(&e, indices);
//                     assert(n_dofs == SizeType(indices.size()));

//                     for (SizeType i = 0; i < n_dofs; ++i) {
//                         idx[i] = indices[i];
//                     }

//                     for (Size_t i = 0; i < nn; ++i) {
//                         for (int di = 0; di < n_var; ++di) {
//                             const Size_t idx_i = i * n_var + di;

//                             for (Size_t j = 0; j < nn; ++j) {
//                                 for (int dj = 0; dj < n_var; ++dj) {
//                                     const Size_t idx_j = j * n_var + dj;
//                                     val[idx_i * n_dofs + idx_j] = element_matrices(elem_idx, idx_i, idx_j);
//                                 }
//                             }
//                         }
//                     }

//                     if (is_block) {
//                         matrix.add_matrix_blocked(idx, idx, val);
//                     } else {
//                         matrix.add_matrix(idx, idx, val);
//                     }
//                 }
//             }
//         }

//         if (mode == SUBTRACT_MODE) {
//             matrix *= -1.0;
//         }

//         UTOPIA_TRACE_REGION_END("LocalToGlobal(libMesh,Intrepid2)");
//     }

//     template <typename Scalar>
//     void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<Scalar>, PetscVector>::apply(
//         const utopia::libmesh::FunctionSpace &space,
//         const LibMeshViewDevice_t<Scalar> &element_vectors,
//         AssemblyMode mode,
//         PetscVector &vector) {
//         assert(false);
//     }

//     template <typename Scalar>
//     void LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<Scalar>, PetscVector>::side_apply(
//         const utopia::libmesh::FunctionSpace &space,
//         const LibMeshViewDevice_t<Scalar> &device_element_vectors,
//         AssemblyMode mode,
//         PetscVector &vector,
//         const std::string &part_name) {
//         UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(libMesh,Intrepid2)");

//         // Type defs
//         using IndexArray_t = Traits<PetscMatrix>::IndexArray;
//         using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
//         using Size_t = Traits<PetscMatrix>::SizeType;
//         using HostView = typename LibMeshViewDevice_t<Scalar>::HostMirror;

//         HostView element_vectors = ::Kokkos::create_mirror_view(device_element_vectors);
//         ::Kokkos::deep_copy(element_vectors, device_element_vectors);

//         if (vector.empty()) {
//             space.create_vector(vector);
//         } else {
//             // Reuse vector
//             if (mode == OVERWRITE_MODE) {
//                 vector *= 0.0;
//             } else if (mode == SUBTRACT_MODE) {
//                 vector *= -1.0;
//             }
//         }

//         auto &m = space.mesh().raw_type();
//         auto &dof_map = space.raw_type_dof_map();

//         {
//             Write<PetscVector> w(vector, utopia::GLOBAL_ADD);

//             auto e_begin = m.active_local_elements_begin();
//             auto e_end = m.active_local_elements_end();

//             if (e_begin != e_end) {
//                 const int n_var = space.n_var();
//                 const SizeType n_dofs = element_vectors.extent(1);
//                 const SizeType nn = n_dofs / n_var;

//                 const bool is_block = vector.block_size() > 1;
//                 IndexArray_t idx(is_block ? nn : n_dofs);
//                 ScalarArray_t val(n_dofs * n_dofs);
//                 std::vector<libMesh::dof_id_type> indices;

//                 int elem_idx = 0;
//                 for (auto e_it = e_begin; e_it != e_end; ++e_it, ++elem_idx) {
//                     const auto &e = **e_it;

//                     dof_map.dof_indices(&e, indices);
//                     assert(n_dofs == SizeType(indices.size()));

//                     for (SizeType i = 0; i < n_dofs; ++i) {
//                         idx[i] = indices[i];
//                     }

//                     for (Size_t i = 0; i < nn; ++i) {
//                         for (int di = 0; di < n_var; ++di) {
//                             const Size_t idx_i = i * n_var + di;
//                             val[idx_i] = element_vectors(elem_idx, idx_i);
//                         }
//                     }

//                     if (is_block) {
//                         vector.add_vector_blocked(idx, val);
//                     } else {
//                         vector.add_vector(idx, val);
//                     }
//                 }
//             }
//         }

//         if (mode == SUBTRACT_MODE) {
//             vector *= -1.0;
//         }

//         UTOPIA_TRACE_REGION_END("LocalToGlobal(libMesh,Intrepid2)");
//     }

//     template class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<double>, PetscMatrix>;
//     template class LocalToGlobal<utopia::libmesh::FunctionSpace, LibMeshViewDevice_t<double>, PetscVector>;

//     template <typename Scalar>
//     void ConvertField<Field<utopia::libmesh::FunctionSpace>, utopia::intrepid2::Field<Scalar>>::apply(
//         const Field<utopia::libmesh::FunctionSpace> &from,
//         utopia::intrepid2::Field<Scalar> &to) {
//         GlobalToLocal<utopia::libmesh::FunctionSpace, Vector, LibMeshViewDevice_t<Scalar>>::apply(
//             *from.space(), from.data(), to.data(), from.tensor_size());

//         to.set_tensor_size(from.tensor_size());
//     }

//     template <typename Scalar>
//     void GlobalToLocal<utopia::libmesh::FunctionSpace,
//                        Traits<utopia::libmesh::FunctionSpace>::Vector,
//                        LibMeshViewDevice_t<Scalar>>::apply(const utopia::libmesh::FunctionSpace &space,
//                                                            const Vector &vector,
//                                                            LibMeshViewDevice_t<Scalar> &device_element_vectors,
//                                                            const int n_comp) {
//         UTOPIA_TRACE_REGION_BEGIN("GlobalToLocal(libMesh,Intrepid2)");

//         using IndexArray_t = Traits<PetscVector>::IndexArray;
//         using ScalarArray_t = Traits<PetscVector>::ScalarArray;

//         assert(vector.comm().size() == 1 || vector.has_ghosts());

//         auto &m = space.mesh().raw_type();
//         auto &dof_map = space.raw_type_dof_map();

//         auto e_begin = m.active_local_elements_begin();
//         auto e_end = m.active_local_elements_end();

//         if (e_begin != e_end) {
//             using HostView = typename LibMeshViewDevice_t<Scalar>::HostMirror;
//             HostView element_vectors = ::Kokkos::create_mirror_view(device_element_vectors);

//             auto num_elem = m.n_active_local_elem();
//             const int n_var = space.n_var();
//             const SizeType n_dofs = element_vectors.extent(1);
//             const SizeType nn = n_dofs / n_var;

//             const bool is_block = vector.block_size() > 1;
//             IndexArray_t idx(is_block ? nn : n_dofs);
//             ScalarArray_t val(n_dofs * n_dofs);
//             std::vector<libMesh::dof_id_type> indices;

//             if (device_element_vectors.extent(0) < num_elem || device_element_vectors.extent(1) < nn) {
//                 device_element_vectors = LibMeshViewDevice_t<Scalar>("Coefficients", num_elem, nn, n_comp);
//             }

//             // FIXME not on device
//             // auto local_view = const_local_view_device(local);

//             Read<PetscVector> r(vector);

//             int elem_idx = 0;
//             for (auto e_it = e_begin; e_it != e_end; ++e_it, ++elem_idx) {
//                 const auto &e = **e_it;

//                 dof_map.dof_indices(&e, indices);
//                 assert(n_dofs == SizeType(indices.size()));

//                 for (SizeType i = 0; i < n_dofs; ++i) {
//                     idx[i] = indices[i];
//                 }

//                 vector.get(idx, val);

//                 for (int i = 0; i < nn; ++i) {
//                     for (int d = 0; d < n_comp; ++d) {
//                         element_vectors(elem_idx, i, d) = val[i * n_comp + d];
//                     }
//                 }
//             }

//             ::Kokkos::deep_copy(device_element_vectors, element_vectors);
//         }

//         UTOPIA_TRACE_REGION_END("GlobalToLocal(libMesh,Intrepid2)");
//     }

//     template class ConvertField<Field<utopia::libmesh::FunctionSpace>, utopia::intrepid2::Field<LibMeshScalar_t>>;
//     template class GlobalToLocal<utopia::libmesh::FunctionSpace,
//                                  Traits<utopia::libmesh::FunctionSpace>::Vector,
//                                  LibMeshViewDevice_t<LibMeshScalar_t>>;

// }  // namespace utopia
