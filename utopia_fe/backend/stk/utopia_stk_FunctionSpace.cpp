#include "utopia_stk_FunctionSpace.hpp"

#include <memory>

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

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return -1; }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return -1; }

    }  // namespace stk
}  // namespace utopia
