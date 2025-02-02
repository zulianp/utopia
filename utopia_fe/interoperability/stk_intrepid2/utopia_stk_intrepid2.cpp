#include "utopia_stk_intrepid2.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc_Matrix.hpp"
#endif

#include "utopia_stk_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_stk_DofMap.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_kokkos_L2Norm.hpp"

// Stk
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Shards_BasicTopologies.hpp>

namespace utopia {

    static ::shards::CellTopology convert_elem_type(::stk::topology::topology_t topo,
                                                    bool no_shell_topologies = false) {
        switch (topo) {
            case ::stk::topology::NODE:
                return ::shards::getCellTopologyData<::shards::Node>();

            case ::stk::topology::LINE_2_1D:
                return ::shards::getCellTopologyData<::shards::Line<>>();
            case ::stk::topology::BEAM_2:
            case ::stk::topology::LINE_2:
                return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Line<>>()
                                            : ::shards::getCellTopologyData<::shards::ShellLine<>>());
            case ::stk::topology::SHELL_LINE_2:
                return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Line<>>()
                                            : ::shards::getCellTopologyData<::shards::ShellLine<>>());
            case ::stk::topology::TRI_3_2D:
                return ::shards::getCellTopologyData<::shards::Triangle<>>();
            case ::stk::topology::TRI_3:
                return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Triangle<>>()
                                            : ::shards::getCellTopologyData<::shards::ShellTriangle<>>());
            case ::stk::topology::SHELL_TRI_3:
                return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Triangle<>>()
                                            : ::shards::getCellTopologyData<::shards::ShellTriangle<>>());
            case ::stk::topology::QUAD_4:
                return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
            case ::stk::topology::QUAD_4_2D:
                return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
            case ::stk::topology::SHELL_QUAD_4:
                return (no_shell_topologies ? ::shards::getCellTopologyData<::shards::Quadrilateral<>>()
                                            : ::shards::getCellTopologyData<::shards::ShellQuadrilateral<>>());
            case ::stk::topology::HEX_8:
                return ::shards::getCellTopologyData<::shards::Hexahedron<>>();
            case ::stk::topology::HEX_27:
                return ::shards::getCellTopologyData<::shards::Hexahedron<27>>();
            case ::stk::topology::TET_4:
                return ::shards::getCellTopologyData<::shards::Tetrahedron<>>();
            default: {
                assert(false);
                return ::shards::getCellTopologyData<::shards::Node>();
            }
        }
    }

    template <typename Scalar>
    class CreateFEFromBuckets {
    public:
        using FE = utopia::intrepid2::FE<Scalar>;
        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;

        static bool apply(const ::stk::mesh::BulkData &bulk_data, const BucketVector_t &buckets, FE &fe, int degree) {
            // assert(buckets.begin() != buckets.end());
            if (buckets.begin() == buckets.end()) {
                // utopia::err() << "[Warning] buckets.begin() == buckets.end()\n";
                return false;
            }

            auto &meta_data = bulk_data.mesh_meta_data();

            auto *first_bucket = *buckets.begin();

            // Dirty hack (FIXME once stk usage is a bit more profficient)
            auto topo = convert_elem_type(first_bucket->topology(), true);
            Size_t n_nodes_x_elem = bulk_data.num_nodes((*first_bucket)[0]);
            Size_t spatial_dim = meta_data.spatial_dimension();
            auto *coords = meta_data.coordinate_field();

            Size_t n_cells = 0;
            for (auto *b_ptr : buckets) {
                n_cells += b_ptr->size();
            }

            StkViewDevice_t<Scalar> device_cell_points("cell_points", n_cells, n_nodes_x_elem, spatial_dim);
            StkIntViewDevice_t device_element_tags("element_tags", n_cells);

            typename StkViewDevice_t<Scalar>::HostMirror cell_points = Kokkos::create_mirror_view(device_cell_points);
            StkIntViewDevice_t::HostMirror element_tags = Kokkos::create_mirror_view(device_element_tags);

            int elem_idx = 0;
            for (const auto &ib : buckets) {
                const Bucket_t &b = *ib;

                const Bucket_t::size_type length = b.size();

                const int block = utopia::stk::extract_set_id_from_bucket(b, b.topology().rank());

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t elem = b[k];
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    element_tags(elem_idx) = block;

                    assert(n_nodes == n_nodes_x_elem);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        const Scalar_t *point = (const Scalar_t *)::stk::mesh::field_data(*coords, node_ids[i]);

                        for (int d = 0; d < spatial_dim; ++d) {
                            cell_points(elem_idx, i, d) = point[d];
                        }
                    }

                    ++elem_idx;
                }
            }

            Kokkos::deep_copy(device_cell_points, cell_points);
            Kokkos::deep_copy(device_element_tags, element_tags);

            fe.init(topo, device_cell_points, degree);
            fe.element_tags() = device_element_tags;
            return true;
        }
    };

    template <typename Scalar>
    void CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::stk::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        using BucketVector_t = ::stk::mesh::BucketVector;

        auto &mesh = space.mesh();
        auto &bulk_data = mesh.bulk_data();
        const BucketVector_t &elem_buckets = utopia::stk::local_elements(bulk_data);

        CreateFEFromBuckets<Scalar>::apply(bulk_data, elem_buckets, fe, degree);
    }

    // Explicit instantation
    template class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar_t>>;

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::stk::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        using BucketVector_t = ::stk::mesh::BucketVector;

        auto &meta_data = space.mesh().meta_data();
        auto &bulk_data = space.mesh().bulk_data();

        // ::stk::mesh::Selector s_universal = meta_data.universal_part();
        ::stk::mesh::Selector selector = meta_data.locally_owned_part();
        const BucketVector_t &side_buckets = bulk_data.get_buckets(meta_data.side_rank(), selector);

        CreateFEFromBuckets<Scalar>::apply(bulk_data, side_buckets, fe, degree);
    }

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::stk::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const std::string &part_name,
        const int degree) {
        using BucketVector_t = ::stk::mesh::BucketVector;

        auto &meta_data = space.mesh().meta_data();
        auto &bulk_data = space.mesh().bulk_data();

        auto part_ptr = meta_data.get_part(part_name);
        if (!part_ptr) {
            part_ptr = meta_data.get_part(stk::SideSet::Cube::convert(part_name));
        }

        assert(part_ptr);
        if (!part_ptr) {
            Utopia::Abort("Unable to find part!");
        }

        ::stk::mesh::Selector part = *part_ptr & meta_data.locally_owned_part();

        const BucketVector_t &side_buckets = bulk_data.get_buckets(meta_data.side_rank(), part);

        CreateFEFromBuckets<Scalar>::apply(bulk_data, side_buckets, fe, degree);
    }

    template class CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<StkScalar_t>>;

#ifdef UTOPIA_ENABLE_PETSC
    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscMatrix>::apply(
        const utopia::stk::FunctionSpace &space,
        const StkViewDevice_t<Scalar> &device_element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix) {
        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Stk,Intrepid2)");

        // Type defs
        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
        using Size_t = Traits<PetscMatrix>::SizeType;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;

        typename StkViewDevice_t<Scalar>::HostMirror element_matrices =
            ::Kokkos::create_mirror_view(device_element_matrices);
        ::Kokkos::deep_copy(element_matrices, device_element_matrices);

        const int n_var = space.n_var();

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Stk,Intrepid2,Matrix handling)");

        if (!matrix.is_block() && n_var != 1) {
            matrix.clear();
        }

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

        // Fix preallocation bug and REMOVE ME
        if (!space.mesh().has_aura() && space.n_var() > 1 && space.comm().size() > 1) {
            MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            m_utopia_warning(
                "using: MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR,"
                "PETSC_FALSE);");
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Stk,Intrepid2,Matrix handling)");

        auto &bulk_data = space.mesh().bulk_data();

        const BucketVector_t &elem_buckets = utopia::stk::local_elements(bulk_data);

        UTOPIA_TRACE_REGION_BEGIN("LocalToGlobal(Stk,Intrepid2,Matrix assembly)");
        {
            Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

            const SizeType n_dofs = element_matrices.extent(1);
            const SizeType nn = n_dofs / n_var;

            IndexArray_t idx(nn);
            ScalarArray_t val(n_dofs * n_dofs);

            const bool is_block = matrix.is_block();

            auto &&local_to_global = space.dof_map().local_to_global();

            Size_t elem_idx = 0;
            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    // get the current node entity and extract the id to fill it into the field
                    Entity_t elem = b[k];
                    // const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;
                    const Size_t n_nodes = bulk_data.num_nodes(elem);
                    UTOPIA_UNUSED(n_nodes);

                    assert(nn == n_nodes);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    if (local_to_global.empty()) {
                        for (Size_t i = 0; i < nn; ++i) {
                            idx[i] = utopia::stk::convert_entity_to_index(node_ids[i]);
                            assert(idx[i] < space.n_dofs());
                            assert(idx[i] >= 0);
                        }
                    } else {
                        for (Size_t i = 0; i < nn; ++i) {
                            idx[i] = local_to_global.block(utopia::stk::convert_entity_to_index(node_ids[i]));
                            assert(idx[i] < space.n_dofs());
                            assert(idx[i] >= 0);
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

                    ++elem_idx;
                }
            }
        }
        UTOPIA_TRACE_REGION_END("LocalToGlobal(Stk,Intrepid2,Matrix assembly)");

        if (mode == SUBTRACT_MODE) {
            matrix *= -1.0;
        }

        UTOPIA_TRACE_REGION_END("LocalToGlobal(Stk,Intrepid2)");
    }

    template <typename Scalar>
    class LocalToGlobalFromBuckets {
    public:
        using FE = utopia::intrepid2::FE<Scalar>;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;

        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;

        static void apply(const utopia::stk::FunctionSpace &space,
                          const BucketVector_t &buckets,
                          const StkViewDevice_t<Scalar> &device_element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector) {
            const int n_var = space.n_var();
            apply(space, buckets, device_element_vectors, mode, vector, n_var);
        }

        static void apply_elemental(const utopia::stk::FunctionSpace &space,
                                    // const BucketVector_t &buckets,
                                    const StkViewDevice_t<Scalar> &device_element_vectors,
                                    AssemblyMode mode,
                                    PetscVector &vector) {
            typename StkViewDevice_t<Scalar>::HostMirror element_vectors =
                ::Kokkos::create_mirror_view(device_element_vectors);
            ::Kokkos::deep_copy(element_vectors, device_element_vectors);

            const Size_t n_cells = element_vectors.extent(0);
            const int tensor_size = element_vectors.extent(1);

            if (empty(vector)) {
                vector.zeros(layout(space.comm(), n_cells * tensor_size, Traits<PetscVector>::determine()));
            } else {
                if (mode == OVERWRITE_MODE) {
                    vector *= 0.0;
                }
            }

            Write<PetscVector> w(vector);

            auto r = range(vector);
            for (Size_t i = 0; i < n_cells; ++i) {
                Size_t offset_i = i * tensor_size;

                for (int t = 0; t < tensor_size; ++t) {
                    vector.add(offset_i + t, element_vectors(i, t));
                }
            }
        }

        static void apply(const utopia::stk::FunctionSpace &space,
                          const BucketVector_t &buckets,
                          const StkViewDevice_t<Scalar> &device_element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector,
                          const int n_var) {
            UTOPIA_TRACE_REGION_BEGIN("LocalToGlobalFromBuckets(Stk,Intrepid2)");

            if (empty(vector)) {
                space.create_vector(vector);
            } else {
                // Ensure the vector has the right block size
                vector.set_block_size(n_var);
                // Reuse vector
                if (mode == OVERWRITE_MODE) {
                    vector *= 0.0;
                } else if (mode == SUBTRACT_MODE) {
                    vector *= -1.0;
                }
            }

            typename StkViewDevice_t<Scalar>::HostMirror element_vectors =
                ::Kokkos::create_mirror_view(device_element_vectors);
            ::Kokkos::deep_copy(element_vectors, device_element_vectors);

            auto &bulk_data = space.mesh().bulk_data();

            {
                Write<PetscVector> w(vector, utopia::GLOBAL_ADD);

                const SizeType n_dofs = element_vectors.extent(1);
                const SizeType nn = n_dofs / n_var;

                IndexArray_t idx(nn);
                ScalarArray_t val(n_dofs);

                const bool is_block = vector.block_size() > 1;

                auto &&local_to_global = space.dof_map().local_to_global();

                Size_t elem_idx = 0;
                for (const auto &ib : buckets) {
                    const Bucket_t &b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        // get the current node entity and extract the id to fill it into the field
                        Entity_t elem = b[k];
                        const Size_t n_nodes = bulk_data.num_nodes(elem);
                        UTOPIA_UNUSED(n_nodes);

                        assert(nn == n_nodes);

                        auto node_ids = bulk_data.begin_nodes(elem);

                        if (local_to_global.empty()) {
                            for (Size_t i = 0; i < nn; ++i) {
                                idx[i] = utopia::stk::convert_entity_to_index(node_ids[i]);
                                assert(idx[i] < space.n_dofs());
                                assert(idx[i] >= 0);
                            }
                        } else {
                            for (Size_t i = 0; i < nn; ++i) {
                                idx[i] = local_to_global.block(utopia::stk::convert_entity_to_index(node_ids[i]));
                                assert(idx[i] < space.n_dofs());
                                assert(idx[i] >= 0);
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

                        ++elem_idx;
                    }
                }
            }

            if (mode == SUBTRACT_MODE) {
                vector *= -1.0;
            }

            // disp(vector);

            UTOPIA_TRACE_REGION_END("LocalToGlobalFromBuckets(Stk,Intrepid2)");
        }
    };

    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::stk::FunctionSpace &space,
        const StkViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector) {
        LocalToGlobalFromBuckets<Scalar>::apply(
            space, utopia::stk::local_elements(space.mesh().bulk_data()), element_vectors, mode, vector);
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscVector>::apply(
        const utopia::stk::FunctionSpace &space,
        const StkViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const int n_var) {
        LocalToGlobalFromBuckets<Scalar>::apply(
            space, utopia::stk::local_elements(space.mesh().bulk_data()), element_vectors, mode, vector, n_var);
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<Scalar>, PetscVector>::side_apply(
        const utopia::stk::FunctionSpace &space,
        const StkViewDevice_t<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &part_name) {
        auto &meta_data = space.mesh().meta_data();
        auto &bulk_data = space.mesh().bulk_data();

        auto part_ptr = meta_data.get_part(part_name);
        if (!part_ptr) {
            part_ptr = meta_data.get_part(stk::SideSet::Cube::convert(part_name));
        }

        assert(part_ptr);
        if (!part_ptr) {
            Utopia::Abort("Unable to find part!");
        }

        ::stk::mesh::Selector part = *part_ptr & meta_data.locally_owned_part();

        LocalToGlobalFromBuckets<Scalar>::apply(
            space, bulk_data.get_buckets(meta_data.side_rank(), part), element_vectors, mode, vector);
    }

    template class LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<StkScalar_t>, PetscMatrix>;
    template class LocalToGlobal<utopia::stk::FunctionSpace, StkViewDevice_t<StkScalar_t>, PetscVector>;

#endif

    template <typename Scalar>
    void ConvertField<Field<utopia::stk::FunctionSpace>, Intrepid2Field<Scalar>>::apply(
        const Field<utopia::stk::FunctionSpace> &from,
        Intrepid2Field<Scalar> &to) {
        GlobalToLocal<utopia::stk::FunctionSpace, Vector, StkViewDevice_t<Scalar>>::apply(
            *from.space(), from.data(), to.data(), from.tensor_size());

        to.set_tensor_size(from.tensor_size());
        to.set_elem_type(from.elem_type());
    }

    template class ConvertField<Field<utopia::stk::FunctionSpace>, Intrepid2Field<StkScalar_t>>;

    template <typename Scalar>
    void ConvertField<Intrepid2Field<Scalar>, Field<utopia::stk::FunctionSpace>>::apply(
        const Intrepid2Field<Scalar> &from,
        Field<utopia::stk::FunctionSpace> &to) {
        if (from.elem_type() == ELEMENT_TYPE) {
            LocalToGlobalFromBuckets<Scalar>::apply_elemental(
                *to.space(),
                // utopia::stk::local_elements(to.space()->mesh().bulk_data()),
                from.data(),
                OVERWRITE_MODE,
                to.data());
        } else {
            assert(false && "IMPLEMENT ME");
            Utopia::Abort("IMPLEMENT ME");
            // LocalToGlobal<utopia::stk::FunctionSpace, Vector, StkViewDevice_t<Scalar>>::apply(
            //     *to.space(), from.data(), OVERWRITE_MODE, to.data(), from.tensor_size());
        }

        to.set_tensor_size(from.tensor_size());
        to.set_elem_type(from.elem_type());
    }

    template class ConvertField<Intrepid2Field<StkScalar_t>, Field<utopia::stk::FunctionSpace>>;

    template <typename Scalar>
    void GlobalToLocal<utopia::stk::FunctionSpace,
                       Traits<utopia::stk::FunctionSpace>::Vector,
                       StkViewDevice_t<Scalar>>::apply(const utopia::stk::FunctionSpace &space,
                                                       const Vector &vector,
                                                       StkViewDevice_t<Scalar> &device_element_vectors,
                                                       const int n_comp) {
        UTOPIA_TRACE_REGION_BEGIN("GlobalToLocal(Stk,Intrepid2)");

        Vector local;
        space.global_to_local(vector, local);

        auto &bulk_data = space.mesh().bulk_data();
        const auto &elem_buckets = utopia::stk::local_elements(bulk_data);
        const Size_t num_elem = utopia::stk::count_local_elements(bulk_data);

        // Dirty hack (FIXME once stk usage is a bit more profficient)
        auto *first_bucket = *elem_buckets.begin();
        auto topo = convert_elem_type(first_bucket->topology(), true);
        Size_t n_nodes_x_elem = bulk_data.num_nodes((*first_bucket)[0]);

        if (Size_t(device_element_vectors.extent(0)) < num_elem ||
            Size_t(device_element_vectors.extent(1)) < n_nodes_x_elem * n_comp) {
            device_element_vectors = StkViewDevice_t<Scalar>("Coefficients", num_elem, n_nodes_x_elem * n_comp);
        }

        typename StkViewDevice_t<Scalar>::HostMirror element_vectors =
            ::Kokkos::create_mirror_view(device_element_vectors);

        // FIXME not on device
        auto local_view = const_local_view_device(local);

        Size_t elem_idx = 0;
        for (const auto &ib : elem_buckets) {
            const auto &b = *ib;
            const SizeType length = b.size();

            for (SizeType k = 0; k < length; ++k) {
                // get the current node entity and extract the id to fill it into the field
                auto elem = b[k];
                const int n_nodes = bulk_data.num_nodes(elem);

                auto node_ids = bulk_data.begin_nodes(elem);

                for (int i = 0; i < n_nodes; ++i) {
                    SizeType node_idx = utopia::stk::convert_entity_to_index(node_ids[i]);
                    for (int d = 0; d < n_comp; ++d) {
                        element_vectors(elem_idx, i * n_comp + d) = local_view.get(node_idx * n_comp + d);
                    }
                }

                ++elem_idx;
            }
        }

        ::Kokkos::deep_copy(device_element_vectors, element_vectors);

        UTOPIA_TRACE_REGION_END("GlobalToLocal(Stk,Intrepid2)");
    }

    template class GlobalToLocal<utopia::stk::FunctionSpace,
                                 Traits<utopia::stk::FunctionSpace>::Vector,
                                 StkViewDevice_t<StkScalar_t>>;

    void l2_norm(const Field<utopia::stk::FunctionSpace> &field, std::vector<StkScalar_t> &norms) {
        using Scalar = typename Traits<utopia::stk::FunctionSpace>::Scalar;
        utopia::kokkos::L2Norm<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>().compute(field, norms);
    }

    // template <class FunctionSpace>
    // auto l2_norm(const utopia::Field<FunctionSpace> &field) -> typename Traits<FunctionSpace>::Scalar {
    //     using Scalar = typename Traits<FunctionSpace>::Scalar;
    //     std::vector<Scalar> results;
    //     L2Norm<FunctionSpace>().compute(field, results);

    //     assert(results.size() == 1);
    //     return results[0];
    // }

    // template <class FunctionSpace>
    // void l2_norm(const utopia::Field<FunctionSpace> &field,
    //              std::vector<typename Traits<FunctionSpace>::Scalar> &results) {

    // }

}  // namespace utopia
