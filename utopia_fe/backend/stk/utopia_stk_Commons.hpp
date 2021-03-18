#ifndef UTOPIA_STK_COMMONS_HPP
#define UTOPIA_STK_COMMONS_HPP

#include "utopia_stk_ForwardDeclarations.hpp"

#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

namespace utopia {
    namespace stk {
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
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_COMMONS_HPP