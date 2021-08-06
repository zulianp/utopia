#include "utopia_petsc_FunctionSpace.hpp"

#include "utopia_petsc_DMDABase.hpp"

namespace utopia {

    namespace petsc {

        class FunctionSpace::Impl : public Configurable {
        public:
            void read(Input &in) override { mesh->read(in); }

            Impl(const Communicator &comm) : mesh(std::make_shared<Mesh>(comm)) {}
            Impl(const std::shared_ptr<Mesh> &mesh) : mesh(mesh) {}
            std::shared_ptr<Mesh> mesh;
        };

        void FunctionSpace::read(Input &in) { impl_->read(in); }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>(mesh)) {}
        FunctionSpace::FunctionSpace(const Communicator &comm) : impl_(utopia::make_unique<Impl>(comm)) {}

        FunctionSpace::~FunctionSpace() = default;

        void FunctionSpace::create_matrix(Matrix &mat) const { impl_->mesh->dm().create_matrix(mat); }
        void FunctionSpace::create_vector(Vector &vec) const { impl_->mesh->dm().create_vector(vec); }

    }  // namespace petsc
}  // namespace utopia
