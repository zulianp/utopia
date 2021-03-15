#include "utopia_moonolith_stk.hpp"

// Utopia fe
#include "utopia_moonolith_Mesh.hpp"
#include "utopia_stk_Mesh.hpp"

// Stk
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

// Std lib
#include <cassert>

namespace utopia {

    void ConvertMesh<utopia::stk::Mesh, utopia::moonolith::Mesh>::apply(const utopia::stk::Mesh &in,
                                                                        utopia::moonolith::Mesh &out) {
        using Bucket_t = ::stk::mesh::Bucket;
        using BucketVector_t = ::stk::mesh::BucketVector;
        using Entity_t = ::stk::mesh::Entity;

        auto &meta_data = in.meta_data();
        auto &bulk_data = in.bulk_data();

        ::stk::mesh::Selector s_universal = meta_data.universal_part();
        const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

        for (const auto &ib : node_buckets) {
            const Bucket_t &b = *ib;
            const Bucket_t::size_type length = b.size();

            for (Bucket_t::size_type k = 0; k < length; ++k) {
                // get the current node entity and extract the id to fill it into the field
                Entity_t node = b[k];
                std::cout << node << std::endl;
            }
        }
    }

}  // namespace utopia
