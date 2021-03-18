#include "utopia_stk_intrepid2.hpp"

#include "utopia_stk_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"
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

    static ::shards::CellTopology convert_elem_type(::stk::topology::topology_t topo) {
        switch (topo) {
            case ::stk::topology::NODE:
                return ::shards::getCellTopologyData<::shards::Node>();
            case ::stk::topology::TRI_3_2D:
                return ::shards::getCellTopologyData<::shards::Triangle<>>();
            case ::stk::topology::SHELL_TRI_3:
                return ::shards::getCellTopologyData<::shards::ShellTriangle<>>();
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
            default:
                return ::shards::getCellTopologyData<::shards::Node>();
        }
    }

    template <typename Scalar>
    void CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>>::apply(
        const utopia::stk::FunctionSpace &space,
        utopia::intrepid2::FE<Scalar> &fe,
        const int degree) {
        // using FunctionSpace = utopia::stk::FunctionSpace;
        using FE = utopia::intrepid2::FE<Scalar>;
        using DynRankView = typename FE::DynRankView;

        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;
        using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
        using Size_t = Traits<utopia::stk::Mesh>::SizeType;

        auto &mesh = space.mesh();
        auto &meta_data = mesh.meta_data();
        auto &bulk_data = mesh.bulk_data();

        ::stk::mesh::Selector s_universal = meta_data.universal_part();
        const BucketVector_t &elem_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, s_universal);

        Size_t n_local_nodes = mesh.n_local_nodes();
        Size_t n_local_elements = mesh.n_local_elements();
        Size_t spatial_dim = mesh.spatial_dimension();
        auto *first_bucket = *elem_buckets.begin();

        // Dirty hack (FIXME once stk usage is a bit more profficient)
        auto topo = convert_elem_type(first_bucket->topology());
        Size_t n_nodes_x_elem = bulk_data.num_nodes((*first_bucket)[0]);
        ::stk::mesh::FieldBase *coords = ::stk::mesh::get_field_by_name("coordinates", meta_data);

        DynRankView cell_points("cell_points", n_local_elements, n_nodes_x_elem, spatial_dim);

        for (const auto &ib : elem_buckets) {
            const Bucket_t &b = *ib;
            const Bucket_t::size_type length = b.size();

            for (Bucket_t::size_type k = 0; k < length; ++k) {
                // get the current node entity and extract the id to fill it into the field
                Entity_t elem = b[k];
                const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;
                const Size_t n_nodes = bulk_data.num_nodes(elem);

                assert(n_nodes == n_nodes_x_elem);

                auto node_ids = bulk_data.begin_nodes(elem);

                for (Size_t i = 0; i < n_nodes; ++i) {
                    const Scalar_t *point = (const Scalar_t *)::stk::mesh::field_data(*coords, node_ids[i]);

                    for (int d = 0; d < spatial_dim; ++d) {
                        cell_points(elem_idx, i, d) = point[d];
                    }
                }
            }
        }

        fe.init(topo, cell_points, degree);
    }

    // Explicit instantations
    template class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>>;

}  // namespace utopia
