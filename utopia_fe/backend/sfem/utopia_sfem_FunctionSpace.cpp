#include "utopia_sfem_FunctionSpace.hpp"

namespace utopia {
    namespace sfem {

        class FunctionSpace::Impl {
        public:
            std::shared_ptr<Mesh> mesh;
            std::string name{"sfem::FunctionSpace"};
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(std::make_unique<Impl>()) {}

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) {
            // TODO
        }

        FunctionSpace::~FunctionSpace() {
            // TODO
        }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) { impl_->mesh = mesh; }

        void FunctionSpace::update(const SimulationTime<Scalar> &) {
            // TODO
        }

        bool FunctionSpace::write(const Path &path, const Vector &x) {
            // TODO
        }

        void FunctionSpace::read(Input &in) {
            if (!impl_->mesh) impl_->mesh = std::make_shared<Mesh>();
            in.require("mesh", *impl_->mesh);
        }

        void FunctionSpace::describe(std::ostream &os) const {
            // TODO
        }

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const {
            // TODO
        }

        const Mesh &FunctionSpace::mesh() const { return *impl_->mesh; }

        Mesh &FunctionSpace::mesh() { return *impl_->mesh; }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const {
            // TODO
        }

        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {
            // TODO
        }

        int FunctionSpace::n_var() const { return 1; }

        void FunctionSpace::create_vector(Vector &v) const {
            // TODO
        }

        void FunctionSpace::create_matrix(Matrix &m) const {
            // TODO
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) const {
            // TODO
        }

        void FunctionSpace::apply_constraints(Vector &v) const {
            // TODO
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) const {
            // TODO
        }

        void FunctionSpace::apply_zero_constraints(Vector &vec) const {
            // TODO
        }

        void FunctionSpace::copy_at_constrained_nodes(const Vector &, Vector &) const {
            // TODO
        }

        void FunctionSpace::overwrite_parts(const std::vector<std::string> &parts,
                                            const std::vector<int> &components,
                                            const Vector &source,
                                            Vector &destination) const {
            // TODO
        }

        void FunctionSpace::set_overwrite_vector(const Vector &v) {
            // TODO
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            // TODO
        }

        bool FunctionSpace::empty() const { return impl_->mesh == nullptr; }

        const std::string &FunctionSpace::name() const { return impl_->name; }
        void FunctionSpace::initialize() {}
        void FunctionSpace::displace(const Vector &displacement) {}
        bool FunctionSpace::read_with_state(Input &in, Field<FunctionSpace> &val) { return false; }
    }  // namespace sfem
}  // namespace utopia
