#include "utopia_petsc_FunctionSpace.hpp"

#include "utopia_FEVar.hpp"
#include "utopia_Options.hpp"
#include "utopia_petsc_DMDABase.hpp"

namespace utopia {

    namespace petsc {

        template <typename T>
        class EmptyUniquePtrDeleter {
        public:
            void operator()(T *) const {}
        };

        class FunctionSpace::Impl : public Configurable {
        public:
            void read(Input &in) override {
                int n_var = 1;
                if (!Options()
                         .add_option("name", name, "Unique name of the function space")
                         .add_option("boundary_conditions", dirichlet_boundary, "Boundary conditions.")
                         .add_option("verbose", verbose, "Verbose output.")
                         .parse(in)) {
                    return;
                }

                variables.read(in);
                mesh->dm().set_default_components(variables.count_variables());
                in.get("mesh", *mesh);

                // auto &&dm = mesh->dm();

                dirichlet_boundary.convert_user_space_names(SideSet::Cube());

                // int i = 0;
                // for (auto &v : variables) {
                //     dm.set_field_name(i++, v.name);
                // }
            }

            Impl(const Communicator &comm) : mesh(std::make_shared<Mesh>(comm)) {}
            Impl(const std::shared_ptr<Mesh> &mesh) : mesh(mesh) {}

            std::shared_ptr<Mesh> mesh;
            DirichletBoundary dirichlet_boundary;
            FEVariables variables;
            std::string name{"main"};
            bool verbose{false};
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
        void FunctionSpace::update(const SimulationTime<Scalar> &time) { impl_->dirichlet_boundary.update(time); }

        std::shared_ptr<FunctionSpace::Mesh> FunctionSpace::mesh_ptr() const {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        int FunctionSpace::n_var() const { return mesh().dm().get_dof(); }

        // Below could be made into a cross-backend interface (given that the same algebra is used)
        bool FunctionSpace::write(const Path &path, const Vector &x) { return mesh().dm().write(path, x); }

        const FunctionSpace::Communicator &FunctionSpace::comm() const {
            // assert(false);
            // Utopia::Abort("IMPLEMENT ME!");
            return mesh().comm();
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
            using IndexSet = Traits::IndexSet;

            auto mesh_view = mesh().view();
            auto r = mesh_view.dof_range();
            int n_var = mesh().dm().get_dof();

            auto &db = impl_->dirichlet_boundary;

            SizeType n_local_dofs = r.extent();
            SizeType n_local_nodes = n_local_dofs / n_var;
            assert(n_local_nodes * n_var == n_local_dofs);

            IndexSet constrains;
            constrains.reserve(n_local_dofs);

            for (auto i = 0; i < n_local_nodes; ++i) {
                for (auto &c : db) {
                    if (mesh_view.is_node_on_boundary(i, c->side)) {
                        const SizeType dof = i * n_var + c->component;
                        constrains.push_back(dof + r.begin());
                    }
                }
            }

            set_zero_rows(m, constrains, 1.);
        }

        void FunctionSpace::apply_constraints(Vector &v) const {
            using IndexSet = Traits::IndexSet;

            auto v_view = local_view_device(v);

            auto mesh_view = mesh().view();
            auto r = mesh_view.dof_range();
            int n_var = mesh().dm().get_dof();

            auto &db = impl_->dirichlet_boundary;

            SizeType n_local_dofs = r.extent();
            SizeType n_local_nodes = n_local_dofs / n_var;
            SizeType local_size = v.local_size();

            assert(local_size == n_local_dofs);
            assert(n_local_nodes * n_var == n_local_dofs);

            for (auto i = 0; i < n_local_nodes; ++i) {
                for (auto &c : db) {
                    if (mesh_view.is_node_on_boundary(i, c->side)) {
                        const SizeType dof = i * n_var + c->component;
                        const Scalar v = c->value();
                        v_view.set(dof, v);
                    }
                }
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) const {
            using IndexSet = Traits::IndexSet;

            auto v_view = local_view_device(v);

            auto mesh_view = mesh().view();
            auto r = mesh_view.dof_range();
            int n_var = mesh().dm().get_dof();

            auto &db = impl_->dirichlet_boundary;

            SizeType n_local_dofs = r.extent();
            SizeType n_local_nodes = n_local_dofs / n_var;
            SizeType local_size = v.local_size();

            assert(local_size == n_local_dofs);
            assert(n_local_nodes * n_var == n_local_dofs);

            IndexSet constrains;
            constrains.reserve(n_local_dofs);

            for (auto i = 0; i < n_local_nodes; ++i) {
                for (auto &c : db) {
                    if (mesh_view.is_node_on_boundary(i, c->side)) {
                        const SizeType dof = i * n_var + c->component;
                        const Scalar v = c->value();
                        v_view.set(dof, v);

                        constrains.push_back(dof + r.begin());
                    }
                }
            }

            set_zero_rows(m, constrains, 1.);
        }

        void FunctionSpace::apply_zero_constraints(Vector &v) const {
            using IndexSet = Traits::IndexSet;

            auto v_view = local_view_device(v);

            auto mesh_view = mesh().view();
            auto r = mesh_view.dof_range();
            int n_var = mesh().dm().get_dof();

            auto &db = impl_->dirichlet_boundary;

            SizeType n_local_dofs = r.extent();
            SizeType n_local_nodes = n_local_dofs / n_var;
            SizeType local_size = v.local_size();

            assert(local_size == n_local_dofs);
            assert(n_local_nodes * n_var == n_local_dofs);

            for (auto i = 0; i < n_local_nodes; ++i) {
                for (auto &c : db) {
                    if (mesh_view.is_node_on_boundary(i, c->side)) {
                        const SizeType dof = i * n_var + c->component;
                        v_view.set(dof, 0.);
                    }
                }
            }
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            // assert(component < n_var());

            DirichletBoundary::UniformCondition dirichlet_boundary{name, value, component};
            dirichlet_boundary.side = SideSet::Cube::convert(name);
            impl_->dirichlet_boundary.add(dirichlet_boundary);
        }

        bool FunctionSpace::empty() const { return mesh().dm().empty(); }

        void FunctionSpace::displace(const Vector &displacement) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        const std::string &FunctionSpace::name() const { return impl_->name; }

        void FunctionSpace::initialize() {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::copy_meta_info_from(const FunctionSpace &other) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void FunctionSpace::create_field(Field<FunctionSpace> &field) {
            auto gv = std::make_shared<Vector>();
            create_vector(*gv);
            field.set_data(gv);
            field.set_space(make_ref(*this));

            field.set_tensor_size(this->n_var());

            // assert(impl_->variables.size() == 1);

            if (!impl_->variables.empty()) {
                field.set_name(impl_->variables[0].name);
                rename(impl_->variables[0].name, *gv);
                field.set_tensor_size(impl_->variables[0].n_components);
            }
        }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {
            mesh().dm().global_to_local(global, local);
        }

        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {
            mesh().dm().local_to_global(local, global);
        }

    }  // namespace petsc
}  // namespace utopia
