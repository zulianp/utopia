#include "utopia_stk_intrepid2.hpp"

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_Matrix.hpp"
#endif

#include "utopia_stk_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_stk_DofMap.hpp"
#include "utopia_stk_FunctionSpace.hpp"

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
            case ::stk::topology::QUAD_4_2D:
                return ::shards::getCellTopologyData<::shards::Quadrilateral<>>();
            case ::stk::topology::SHELL_QUAD_4:
                return ::shards::getCellTopologyData<::shards::ShellQuadrilateral<>>();
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
        using DynRankView = typename FE::DynRankView;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;

        static bool apply(const ::stk::mesh::BulkData &bulk_data, const BucketVector_t &buckets, FE &fe, int degree) {
            assert(buckets.begin() != buckets.end());
            if (buckets.begin() == buckets.end()) return false;

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

            DynRankView cell_points("cell_points", n_cells, n_nodes_x_elem, spatial_dim);

            int elem_idx = 0;
            for (const auto &ib : buckets) {
                const Bucket_t &b = *ib;

                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    Entity_t elem = b[k];
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

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

            fe.init(topo, cell_points, degree);
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
    template class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>>;

    template <typename Scalar>
    void CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::stk::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        using BucketVector_t = ::stk::mesh::BucketVector;

        auto &meta_data = space.mesh().meta_data();
        auto &bulk_data = space.mesh().bulk_data();

        ::stk::mesh::Selector s_universal = meta_data.universal_part();
        const BucketVector_t &side_buckets = bulk_data.get_buckets(meta_data.side_rank(), s_universal);

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

        ::stk::mesh::Selector s_universal = *meta_data.get_part(part_name);
        const BucketVector_t &side_buckets = bulk_data.get_buckets(meta_data.side_rank(), s_universal);

        CreateFEFromBuckets<Scalar>::apply(bulk_data, side_buckets, fe, degree);
    }

    template class CreateFEOnBoundary<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>>;

#ifdef UTOPIA_WITH_PETSC
    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscMatrix>::apply(
        const utopia::stk::FunctionSpace &space,
        const ::Kokkos::DynRankView<Scalar> &element_matrices,
        AssemblyMode mode,
        PetscMatrix &matrix) {
        // Type defs
        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;
        using Size_t = Traits<PetscMatrix>::SizeType;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;

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
        if (space.n_var() > 1 && space.comm().size() > 1) {
            MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            m_utopia_warning(
                "using: MatSetOption(matrix.raw_type(), MAT_NEW_NONZERO_ALLOCATION_ERR,"
                "PETSC_FALSE);");
        }

        auto &bulk_data = space.mesh().bulk_data();

        const BucketVector_t &elem_buckets = utopia::stk::local_elements(bulk_data);

        {
            Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

            const int n_var = space.n_var();
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

        if (mode == SUBTRACT_MODE) {
            matrix *= -1.0;
        }
    }

    template <typename Scalar>
    class LocalToGlobalFromBuckets {
    public:
        using FE = utopia::intrepid2::FE<Scalar>;
        using DynRankView = typename FE::DynRankView;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;

        using IndexArray_t = Traits<PetscMatrix>::IndexArray;
        using ScalarArray_t = Traits<PetscMatrix>::ScalarArray;

        static void apply(const utopia::stk::FunctionSpace &space,
                          const BucketVector_t &buckets,
                          const ::Kokkos::DynRankView<Scalar> &element_vectors,
                          AssemblyMode mode,
                          PetscVector &vector) {
            //

            if (empty(vector)) {
                space.create_vector(vector);
            } else {
                // Reuse vector
                if (mode == OVERWRITE_MODE) {
                    vector *= 0.0;
                } else if (mode == SUBTRACT_MODE) {
                    vector *= -1.0;
                }
            }

            auto &bulk_data = space.mesh().bulk_data();

            {
                Write<PetscVector> w(vector, utopia::GLOBAL_ADD);

                const int n_var = space.n_var();
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
        }
    };

    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscVector>::apply(
        const utopia::stk::FunctionSpace &space,
        const ::Kokkos::DynRankView<Scalar> &element_vectors,
        AssemblyMode mode,
        PetscVector &vector) {
        LocalToGlobalFromBuckets<Scalar>::apply(
            space, utopia::stk::local_elements(space.mesh().bulk_data()), element_vectors, mode, vector);
    }

    template <typename Scalar>
    void LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<Scalar>, PetscVector>::side_apply(
        const utopia::stk::FunctionSpace &space,
        const DynRankView &element_vectors,
        AssemblyMode mode,
        PetscVector &vector,
        const std::string &part_name) {
        auto &meta_data = space.mesh().meta_data();
        auto &bulk_data = space.mesh().bulk_data();

        ::stk::mesh::Selector s_universal = *meta_data.get_part(part_name);

        LocalToGlobalFromBuckets<Scalar>::apply(
            space, bulk_data.get_buckets(meta_data.side_rank(), s_universal), element_vectors, mode, vector);
    }

    template class LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<double>, PetscMatrix>;
    template class LocalToGlobal<utopia::stk::FunctionSpace, ::Kokkos::DynRankView<double>, PetscVector>;

#endif

}  // namespace utopia
