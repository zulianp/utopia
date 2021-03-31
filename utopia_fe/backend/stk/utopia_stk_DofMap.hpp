#ifndef UTOPIA_STK_DOF_MAP_HPP
#define UTOPIA_STK_DOF_MAP_HPP

#include "utopia_Describable.hpp"

#include "utopia_Traits.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

namespace utopia {
    namespace stk {

        class DofMap final : public Describable {
        public:
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Communicator = Traits<FunctionSpace>::Communicator;

            void init(::stk::mesh::BulkData &bulk_data);
            DofMap();
            ~DofMap();

            void describe(std::ostream &os) const override;

            bool empty() const;

            const IndexArray &local_to_global() const;
            const IndexArray &d_nnz() const;
            const IndexArray &o_nnz() const;

        private:
            class Impl;
            class DofExchange;
            std::unique_ptr<Impl> impl_;

            void init_parallel(const Communicator &comm, ::stk::mesh::BulkData &bulk_data);
            void init_serial(::stk::mesh::BulkData &bulk_data);
        };
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_DOF_MAP_HPP
