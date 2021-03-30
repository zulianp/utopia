#include "utopia_stk_Commons.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace utopia {
    namespace stk {

        size_t count_entities(const ::stk::mesh::BulkData &bulk_data,
                              const ::stk::topology::rank_t &rank,
                              const ::stk::mesh::Selector &selector) {
            const auto &buckets = bulk_data.get_buckets(rank, selector);
            size_t ret = 0;
            std::vector<int> procs;
            for (auto *b_ptr : buckets) {
                const auto &b = *b_ptr;
                ret += b.size();
            }

            return ret;
        }

        size_t count_local_nodes(const ::stk::mesh::BulkData &bulk_data) {
            return count_entities(
                bulk_data, ::stk::topology::NODE_RANK, bulk_data.mesh_meta_data().locally_owned_part());
        }
        size_t count_local_elements(const ::stk::mesh::BulkData &bulk_data) {
            return count_entities(
                bulk_data, ::stk::topology::ELEMENT_RANK, bulk_data.mesh_meta_data().locally_owned_part());
        }

        size_t count_universal_nodes(const ::stk::mesh::BulkData &bulk_data) {
            return count_entities(bulk_data, ::stk::topology::NODE_RANK, bulk_data.mesh_meta_data().universal_part());
        }
        size_t count_universal_elements(const ::stk::mesh::BulkData &bulk_data) {
            return count_entities(
                bulk_data, ::stk::topology::ELEMENT_RANK, bulk_data.mesh_meta_data().universal_part());
        }

        const ::stk::mesh::BucketVector &shared_nodes(const ::stk::mesh::BulkData &bulk_data) {
            return bulk_data.get_buckets(::stk::topology::NODE_RANK, bulk_data.mesh_meta_data().globally_shared_part());
        }

        const ::stk::mesh::BucketVector &local_nodes(const ::stk::mesh::BulkData &bulk_data) {
            return bulk_data.get_buckets(::stk::topology::NODE_RANK, bulk_data.mesh_meta_data().locally_owned_part());
        }

        const ::stk::mesh::BucketVector &local_elements(const ::stk::mesh::BulkData &bulk_data) {
            return bulk_data.get_buckets(::stk::topology::ELEMENT_RANK,
                                         bulk_data.mesh_meta_data().locally_owned_part());
        }

        const ::stk::mesh::BucketVector &universal_nodes(const ::stk::mesh::BulkData &bulk_data) {
            return bulk_data.get_buckets(::stk::topology::NODE_RANK, bulk_data.mesh_meta_data().universal_part());
        }
        const ::stk::mesh::BucketVector &universal_elements(const ::stk::mesh::BulkData &bulk_data) {
            return bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, bulk_data.mesh_meta_data().universal_part());
        }

    }  // namespace stk
}  // namespace utopia