#ifndef UTOPIA_STK_DOF_MAP_HPP
#define UTOPIA_STK_DOF_MAP_HPP

#include "utopia_Describable.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"

#include <memory>

namespace utopia {
    namespace stk {

        class DofMap final : public Describable {
        public:
            void init(::stk::mesh::BulkData &bulk_data);
            DofMap();
            ~DofMap();

            void describe(std::ostream &os) const override;

            bool empty() const;

        private:
            class Impl;
            class DofExchange;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_DOF_MAP_HPP
