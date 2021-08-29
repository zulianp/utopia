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

        const FunctionSpace::Mesh &FunctionSpace::mesh() const {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        FunctionSpace::Mesh &FunctionSpace::mesh() {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>(mesh)) {}
        FunctionSpace::FunctionSpace(const Communicator &comm) : impl_(utopia::make_unique<Impl>(comm)) {}

        FunctionSpace::~FunctionSpace() = default;

        void FunctionSpace::create_matrix(Matrix &mat) const { impl_->mesh->dm().create_matrix(mat); }
        void FunctionSpace::create_vector(Vector &vec) const { impl_->mesh->dm().create_vector(vec); }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) {}

        std::shared_ptr<FunctionSpace::Mesh> FunctionSpace::mesh_ptr() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        int FunctionSpace::n_var() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        // Below could be made into a cross-backend interface (given that the same algebra is used)
        bool FunctionSpace::write(const Path &path, const Vector &x) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        const FunctionSpace::Communicator &FunctionSpace::comm() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::apply_constraints(Vector &v) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::apply_zero_constraints(Vector &vec) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        bool FunctionSpace::empty() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::displace(const Vector &displacement) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        const std::string &FunctionSpace::name() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::initialize() {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::create_field(Field<FunctionSpace> &field) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

    }  // namespace petsc
}  // namespace utopia
