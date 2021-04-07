
#ifndef UTOPIA_STK_COMMONS_HPP
#define UTOPIA_STK_COMMONS_HPP

#include "utopia_stk_ForwardDeclarations.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_topology/topology.hpp>

namespace utopia {
    namespace stk {

        template <typename Int>
        inline auto convert_stk_index_to_index(Int idx) -> Int {
            return idx - 1;
        }

        template <typename Int>
        inline auto convert_index_to_stk_index(Int idx) -> Int {
            return idx + 1;
        }

        inline auto convert_entity_to_index(const ::stk::mesh::Entity &entity) { return entity.local_offset() - 1; }

        inline int manifold_dim(::stk::topology::topology_t topo) {
            switch (topo) {
                case ::stk::topology::NODE:
                    return 0;
                case ::stk::topology::TRI_3_2D:
                case ::stk::topology::SHELL_TRI_3:
                case ::stk::topology::QUAD_4_2D:
                case ::stk::topology::SHELL_QUAD_4:
                    return 2;
                case ::stk::topology::HEX_8:
                case ::stk::topology::HEX_27:
                case ::stk::topology::TET_4:
                    return 3;
                default:
                    return -1;
            }
        }

        int extract_sideset_id(const std::string &name);

        size_t count_entities(const ::stk::mesh::BulkData &bulk_data,
                              const ::stk::topology::rank_t &rank,
                              const ::stk::mesh::Selector &selector);

        size_t count_local_nodes(const ::stk::mesh::BulkData &bulk_data);
        size_t count_local_elements(const ::stk::mesh::BulkData &bulk_data);

        size_t count_universal_nodes(const ::stk::mesh::BulkData &bulk_data);
        size_t count_universal_elements(const ::stk::mesh::BulkData &bulk_data);

        const ::stk::mesh::BucketVector &local_nodes(const ::stk::mesh::BulkData &bulk_data);
        const ::stk::mesh::BucketVector &local_elements(const ::stk::mesh::BulkData &bulk_data);

        const ::stk::mesh::BucketVector &universal_nodes(const ::stk::mesh::BulkData &bulk_data);
        const ::stk::mesh::BucketVector &universal_elements(const ::stk::mesh::BulkData &bulk_data);

        const ::stk::mesh::BucketVector &shared_nodes(const ::stk::mesh::BulkData &bulk_data);

    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_COMMONS_HPP