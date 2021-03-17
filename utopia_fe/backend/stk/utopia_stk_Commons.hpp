#ifndef UTOPIA_STK_COMMONS_HPP
#define UTOPIA_STK_COMMONS_HPP

#include "utopia_stk_ForwardDeclarations.hpp"

#include <stk_mesh/base/Entity.hpp>

namespace utopia {
    namespace stk {
        inline static auto convert_entity_to_index(const ::stk::mesh::Entity &entity) {
            return entity.local_offset() - 1;
        }
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_COMMONS_HPP