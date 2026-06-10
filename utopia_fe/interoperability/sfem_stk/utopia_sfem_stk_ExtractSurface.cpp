#include "utopia_sfem_stk_ExtractSurface.hpp"

#include "utopia_sfem_Mesh.hpp"

// Utopia::Stk
#include "utopia_stk_Commons.hpp"
#include "utopia_stk_Mesh.hpp"

// Stk
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Entity.hpp>
#if TRILINOS_MAJOR_MINOR_VERSION < 160200
#include <stk_mesh/base/GetBuckets.hpp>
#endif
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "smesh_mesh.hpp"

using BucketVector_t = ::stk::mesh::BucketVector;
using Bucket_t = ::stk::mesh::Bucket;
using Entity_t = ::stk::mesh::Entity;
using Scalar_t = utopia::Traits<utopia::stk::Mesh>::Scalar;
using Size_t = utopia::Traits<utopia::stk::Mesh>::SizeType;
using MetaData_t = ::stk::mesh::MetaData;
using BulkData_t = ::stk::mesh::BulkData;

namespace utopia {
    using idx_t = ::smesh::idx_t;
    using geom_t = ::smesh::geom_t;

    static enum ::smesh::ElemType convert_elem_type(::stk::topology::topology_t topo) {
        switch (topo) {
            case ::stk::topology::NODE:
                return ::smesh::NODE1;
            case ::stk::topology::LINE_2:
            case ::stk::topology::BEAM_2:
                return ::smesh::EDGE2;
            case ::stk::topology::TRI_3:
            case ::stk::topology::TRI_3_2D:
                return ::smesh::TRI3;
            case ::stk::topology::SHELL_TRI_3:
                return ::smesh::TRISHELL3;
            case ::stk::topology::QUAD_4:
            case ::stk::topology::QUAD_4_2D:
            case ::stk::topology::SHELL_QUAD_4:
                return ::smesh::QUAD4;
            case ::stk::topology::HEX_8:
                return ::smesh::HEX8;
            // case ::stk::topology::HEX_27:
            //     return HEX27;
            case ::stk::topology::TET_4:
                return ::smesh::TET4;
            default: {
                assert(false);
                Utopia::Abort("Element type not supported!");
                return ::smesh::INVALID;
            }
        }
    }

    static void extract_selection(const ::stk::mesh::BulkData &bulk_data,
                                  const ::stk::mesh::Selector &selector,
                                  const ::stk::topology::rank_t &topo,
                                  utopia::sfem::Mesh &out) {
        auto &meta_data = bulk_data.mesh_meta_data();
        auto *coords = meta_data.coordinate_field();
        const int spatial_dim = meta_data.spatial_dimension();

        const bool has_aura = bulk_data.is_automatic_aura_on();

        const BucketVector_t &elem_buckets = bulk_data.get_buckets(topo, selector);

        Bucket_t::size_type n_selected_elements = 0;
        Bucket_t::size_type n_selected_nodes = 0;

        Bucket_t::size_type n_universal_nodes = utopia::stk::count_universal_nodes(bulk_data);

        std::vector<SizeType> node_mapping(n_universal_nodes, -1);

        int nodesxelem = 0;

        for (const auto &b_ptr : elem_buckets) {
            auto &b = *b_ptr;

            const Bucket_t::size_type length = b.size();

            n_selected_elements += length;

            for (Bucket_t::size_type k = 0; k < length; ++k) {
                Entity_t elem = b[k];

                const Size_t n_nodes = bulk_data.num_nodes(elem);
                auto node_ids = bulk_data.begin_nodes(elem);

                nodesxelem = n_nodes;

                for (Size_t i = 0; i < n_nodes; ++i) {
                    auto idx = utopia::stk::convert_entity_to_index(node_ids[i]);

                    if (node_mapping[idx] == -1) {
                        node_mapping[idx] = n_selected_nodes++;
                    }
                }
            }
        }

        // std::cout << "Num elem: " << n_selected_elements << "\n";

        ////////////////////////////////////////////////////////////////
        auto elements_buffer = ::smesh::create_host_buffer<idx_t>(nodesxelem, n_selected_elements);
        auto points_buffer = ::smesh::create_host_buffer<geom_t>(spatial_dim, n_selected_nodes);
        auto node_mapping_buffer = ::smesh::create_host_buffer<idx_t>(n_selected_nodes);

        auto elements = elements_buffer->data();
        auto out_points = points_buffer->data();
        auto out_node_mapping = node_mapping_buffer->data();

        ////////////////////////////////////////////////////////////////
        // Nodes
        ////////////////////////////////////////////////////////////////

        const BucketVector_t &node_buckets =
            bulk_data.get_buckets(::stk::topology::NODE_RANK, meta_data.universal_part());

        for (const auto &ib : node_buckets) {
            const Bucket_t &b = *ib;
            const Bucket_t::size_type length = b.size();

            for (Bucket_t::size_type k = 0; k < length; ++k) {
                Entity_t node = b[k];
                if (has_aura && bulk_data.in_receive_ghost(node)) continue;

                auto reduced_index = node_mapping[utopia::stk::convert_entity_to_index(node)];

                if (reduced_index == -1) continue;

                assert(Bucket_t::size_type(reduced_index) < n_selected_nodes);

                const Scalar_t *points = (const Scalar_t *)::stk::mesh::field_data(*coords, node);

                for (int d = 0; d < spatial_dim; ++d) {
                    out_points[d][reduced_index] = points[d];
                }

                out_node_mapping[reduced_index] = utopia::stk::convert_entity_to_index(node);
            }
        }

        ////////////////////////////////////////////////////////////////
        // Elements
        ////////////////////////////////////////////////////////////////

        SizeType selected_elem_idx = 0;
        enum ::smesh::ElemType element_type = ::smesh::INVALID;
        for (const auto &ib : elem_buckets) {
            const Bucket_t &b = *ib;

            auto sfem_type = convert_elem_type(b.topology());
            element_type = sfem_type;
            const Bucket_t::size_type length = b.size();

            for (Bucket_t::size_type k = 0; k < length; ++k) {
                Entity_t elem = b[k];
                const Size_t n_nodes = bulk_data.num_nodes(elem);
                auto node_ids = bulk_data.begin_nodes(elem);

                for (Size_t i = 0; i < n_nodes; ++i) {
                    elements[i][selected_elem_idx] = node_mapping[utopia::stk::convert_entity_to_index(node_ids[i])];
                }

                selected_elem_idx++;
            }
        }

        auto mesh = std::make_shared<::smesh::Mesh>(
            ::smesh::Communicator::wrap(out.comm().get()), element_type, elements_buffer, points_buffer);
        mesh->set_node_mapping(node_mapping_buffer);
        out.set_raw_type(mesh);
    }

    // void ExtractSurface<utopia::stk::Mesh, utopia::sfem::Mesh>::apply(const utopia::stk::Mesh &in,
    //                                                                   utopia::sfem::Mesh &out) {
    //     UTOPIA_TRACE_SCOPE("ExtractSurface::apply");

    //     auto &meta_data = in.meta_data();
    //     auto &bulk_data = in.bulk_data();
    //     auto topo = meta_data.side_rank();

    //     extract_selection(bulk_data, meta_data.locally_owned_part(), topo, out);
    // }

    void ExtractSurface<utopia::stk::Mesh, utopia::sfem::Mesh>::apply(const utopia::stk::Mesh &in,
                                                                      utopia::sfem::Mesh &out,
                                                                      const std::vector<std::string> &exclude) {
        UTOPIA_TRACE_SCOPE("ExtractSurface::apply");

        auto &meta_data = in.meta_data();
        auto &bulk_data = in.bulk_data();
        auto topo = meta_data.side_rank();

        if (exclude.empty()) {
            extract_selection(bulk_data, meta_data.locally_owned_part(), topo, out);
            return;

        } else {
            ::stk::mesh::Selector selector = meta_data.locally_owned_part();

            for (auto &name : exclude) {
                auto *part = meta_data.get_part(name);
                if (part) {
                    selector &= !*part;
                }
            }

            extract_selection(bulk_data, selector, topo, out);
        }
    }

}  // namespace utopia
