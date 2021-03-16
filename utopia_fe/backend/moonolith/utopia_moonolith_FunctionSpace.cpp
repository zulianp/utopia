#include "utopia_moonolith_FunctionSpace.hpp"

#include <memory>

namespace utopia {
    namespace moonolith {

        template <int Dim>
        using MoonolithFunctionSpace = ::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, Dim>>;

        class FunctionSpace::Impl {
        public:
            template <int Dim>
            void wrap(const std::shared_ptr<MoonolithFunctionSpace<Dim>> &ptr) {
                space = ptr;
            }

            std::shared_ptr<Mesh> mesh;
            std::shared_ptr<void> space;
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
        }

        FunctionSpace::~FunctionSpace() = default;

        void FunctionSpace::write(const Path &path, const Vector &x) {}
        void FunctionSpace::read(Input &in) {}
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

        template <int Dim>
        ::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, Dim>> &FunctionSpace::raw_type() {}

        FunctionSpace::SizeType FunctionSpace::n_dofs() const {}
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {}

        // Explicit instantiations
        template ::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 1>> &FunctionSpace::raw_type();
        template ::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 2>> &FunctionSpace::raw_type();
        template ::moonolith::FunctionSpace<::moonolith::Mesh<FunctionSpace::Scalar, 3>> &FunctionSpace::raw_type();

    }  // namespace moonolith
}  // namespace utopia
