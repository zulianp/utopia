#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

// All stk includes
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace utopia {
    namespace stk {

        class FunctionSpace::Impl {
        public:
            std::shared_ptr<Mesh> mesh;
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

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->mesh->n_nodes(); }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return impl_->mesh->n_local_nodes(); }

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

            for (auto *b_ptr : node_buckets) {
                auto &b = *b_ptr;
            }
        }

    }  // namespace stk
}  // namespace utopia
