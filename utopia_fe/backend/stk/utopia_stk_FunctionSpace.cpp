#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_Options.hpp"

#include "utopia_stk_Commons.hpp"

// All stk includes
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <memory>
#include <unordered_set>

namespace utopia {
    namespace stk {

        class FunctionSpace::Impl {
        public:
            std::shared_ptr<Mesh> mesh;
            int n_var{1};
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
        }

        FunctionSpace::~FunctionSpace() = default;

        bool FunctionSpace::write(const Path &path, const Vector &x) { return false; }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) { impl_->mesh = mesh; }

        void FunctionSpace::read(Input &in) {
            auto mesh = std::make_shared<Mesh>();
            mesh->read(in);
            init(mesh);

            if (!Options().add_option("n_var", impl_->n_var, "Number of variables per node.").parse(in)) {
                return;
            }
        }

        void FunctionSpace::describe(std::ostream &os) const {}

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }

        const Mesh &FunctionSpace::mesh() const {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        Mesh &FunctionSpace::mesh() {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->mesh->n_nodes() * impl_->n_var; }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {
            return impl_->mesh->n_local_nodes() * impl_->n_var;
        }

        void FunctionSpace::create_vector(Vector &v) { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }

        void FunctionSpace::create_local_vector(Vector &v) {
            // FIXME
            v.zeros(layout(comm(), n_local_dofs(), n_dofs()));
        }

        void FunctionSpace::create_matrix(Matrix &m) {
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;

            auto vl = layout(comm(), n_local_dofs(), n_dofs());
            auto ml = square_matrix_layout(vl);

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            const ::stk::mesh::Selector selector = meta_data.locally_owned_part();
            const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, selector);

            // Slow but for the moment it will do...
            const SizeType nln = mesh().n_local_nodes();
            std::vector<std::unordered_set<SizeType>> node2node(nln);

            for (auto *b_ptr : node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    const Entity_t elem = b[k];
                    const auto node_ids = bulk_data.begin_nodes(elem);
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        const auto node_i = utopia::stk::convert_entity_to_index(node_ids[i]);

                        for (Size_t j = 0; j < n_nodes; ++j) {
                            const auto node_j = utopia::stk::convert_entity_to_index(node_ids[j]);
                            node2node[node_i].insert(node_j);
                        }
                    }
                }
            }

            IndexSet d_nnz(nln, 0), o_nnz(nln, 0);
            for (SizeType i = 0; i < nln; ++i) {
                d_nnz[i] = node2node[i].size();
            }

            if (impl_->n_var == 1) {
                m.sparse(ml, d_nnz, o_nnz);
            } else {
                m.block_sparse(ml, d_nnz, o_nnz, impl_->n_var);
            }
        }

    }  // namespace stk
}  // namespace utopia
