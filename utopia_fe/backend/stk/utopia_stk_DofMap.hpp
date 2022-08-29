#ifndef UTOPIA_STK_DOF_MAP_HPP
#define UTOPIA_STK_DOF_MAP_HPP

#include "utopia_Describable.hpp"

#include "utopia_Traits.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

#include "utopia_GlobalIndex.hpp"
#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

namespace utopia {

    template <>
    class GlobalIndex<utopia::stk::FunctionSpace, DYNAMIC_SIZE> {
    public:
        using IndexArray = Traits<utopia::stk::FunctionSpace>::IndexArray;
        using SizeType = Traits<utopia::stk::FunctionSpace>::SizeType;

        GlobalIndex(const IndexArray &idx, const int n_var = 1) : idx_(idx), n_var_(n_var) {}

        inline SizeType operator()(const SizeType local_idx) const {
            assert(n_var_ == 1);
            return idx_[local_idx];
        }

        inline SizeType global_node_id(const SizeType local_idx) const { return idx_[local_idx]; }

        inline SizeType operator()(const SizeType local_idx, const int component) const {
            assert(component < n_var_);
            assert(local_idx < static_cast<SizeType>(idx_.size()));
            return idx_[local_idx] * n_var_ + component;
        }

        inline SizeType block(const SizeType local_idx) {
            assert(local_idx < static_cast<SizeType>(idx_.size()));
            return idx_[local_idx];
        }

        inline bool empty() const { return idx_.empty(); }
        inline auto size() const { return idx_.size(); }

        const IndexArray &idx_;
        int n_var_;
    };

    namespace stk {

        class DofMap final : public Describable {
        public:
            using IndexArray = Traits<FunctionSpace>::IndexArray;
            using Mesh = Traits<FunctionSpace>::Mesh;
            using SizeType = Traits<FunctionSpace>::SizeType;
            using Communicator = Traits<FunctionSpace>::Communicator;
            using Vector = Traits<FunctionSpace>::Vector;
            using GlobalIndex = utopia::GlobalIndex<FunctionSpace>;

            void set_valid_local_id_mode(const bool val);
            void init(Mesh &mesh, const bool verbose = false);

            DofMap();
            ~DofMap();

            void describe(std::ostream &os) const override;

            bool empty() const;

            GlobalIndex local_to_global() const;
            const IndexArray &d_nnz() const;
            const IndexArray &o_nnz() const;
            SizeType max_nnz_x_row() const;

            int n_var() const;
            void set_n_var(const int n_var);

            void global_to_local(const Vector &global, Vector &local) const;
            void local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const;

            SizeType shift_aura_idx(const SizeType idx) const;

        private:
            class Impl;
            class DofExchange;
            std::unique_ptr<Impl> impl_;

            void init_parallel(const Communicator &comm, ::stk::mesh::BulkData &bulk_data);
            void init_serial(::stk::mesh::BulkData &bulk_data);

            void init(::stk::mesh::BulkData &bulk_data);

            void exchange_shared_dofs(const Communicator &comm,
                                      ::stk::mesh::BulkData &bulk_data,
                                      bool add_to_nnz_pattern = true);

            void describe_mesh_connectivity(Mesh &mesh);
            void describe_debug(Mesh &mesh, std::ostream &os) const;
        };
    }  // namespace stk
}  // namespace utopia

#endif  // UTOPIA_STK_DOF_MAP_HPP
