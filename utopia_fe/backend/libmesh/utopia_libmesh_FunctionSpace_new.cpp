#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_Options.hpp"

// All libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/reference_counter.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

namespace utopia {
    namespace libmesh {

        class FunctionSpace::Impl {
        public:
            std::shared_ptr<Mesh> mesh;
            std::shared_ptr<libMesh::EquationSystems> systems;
        };

        class FunctionSpace::Var : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("fe_family", fe_family, "The finite element family from enum libMesh::FEFamily.")
                         .add_option("order", order, "The finite element order from enum libMesh::Order.")
                         .add_option("var", name, "Name of the variable")
                         .add_option("components", components, "number of vector components")
                         .parse(in)) {
                    return;
                }
            }

            void add_to_system(libMesh::System &sys) {
                if (components == 1) {
                    sys.add_variable(name,
                                     libMesh::Utility::string_to_enum<libMesh::Order>(order),
                                     libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family));
                } else {
                    const int n = std::min(components, N_SUFFIXES);
                    for (int i = 0; i < n; ++i) {
                        sys.add_variable(name + suffixes[i],
                                         libMesh::Utility::string_to_enum<libMesh::Order>(order),
                                         libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family));
                    }

                    for (int i = N_SUFFIXES; i < components; ++i) {
                        sys.add_variable(name + "_" + std::to_string(i),
                                         libMesh::Utility::string_to_enum<libMesh::Order>(order),
                                         libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family));
                    }
                }
            }

            void describe(std::ostream &os) const override { os << "Var"; }

            std::string fe_family{"LAGRANGE"};
            std::string order{"FIRST"};
            std::string name{"var"};
            int components{1};

            static constexpr const int N_SUFFIXES = 4;
            static constexpr const char *suffixes[N_SUFFIXES] = {"_x", "_y", "_z", "_t"};
        };

        class FunctionSpace::BC : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                std::string bc_name = "bc";
                std::string bc_type = "dirichlet";
                int var = 0;
                std::string type = "constant";
                std::string value = "0";

                if (!Options()
                         .add_option("bc_name", bc_name, "Name of boundary condition.")
                         .add_option("bc_type", bc_type, "Type of boundary condition.")
                         .add_option("var", var, "Variable number associated with boundary condition.")
                         .add_option("type", type, "Type of function: constant|expr.")
                         .add_option("value", value, "Value of function")
                         .parse(in)) {
                    return;
                }
            }

            void add_to_system(libMesh::System &sys) {
                // std::vector<unsigned int> vars(1);
                // vars[0] = cond.expr().left().space_ptr()->subspace_id();
                // std::set<libMesh::boundary_id_type> bt;
                // bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
                // libMesh::DirichletBoundary d_bc(
                //     bt, vars, LibMeshLambdaFunction<libMesh::Real>(cond.expr().right().fun()));
                // cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
            }

            void describe(std::ostream &os) const override { os << "BC"; }
        };

        class FunctionSpace::Sys : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("system_name", system_name, "Name of system.")
                         .add_option("system_type",
                                     system_type,
                                     "Type of system associated with function space: linear_implicit|explicit.")
                         .parse(in)) {
                    return;
                }

                in.get("variables", [&](Input &in) {
                    in.get_all([&](Input &in) {
                        Var v;
                        v.read(in);
                        vars.push_back(std::move(v));
                    });
                });

                in.get("boundary_conditions", [&](Input &in) {
                    in.get_all([&](Input &in) {
                        BC bc;
                        bc.read(in);
                        bcs.push_back(std::move(bc));
                    });
                });

                if (vars.empty()) {
                    // Users did not care, so lets give them a scalar linear fe
                    vars.push_back(Var());
                }
            }

            void add_to_space(Impl &impl) {
                if ("linear_implicit" == system_type) {
                    impl.systems->add_system<libMesh::LinearImplicitSystem>(system_name);
                } else if ("explicit" == system_type) {
                    impl.systems->add_system<libMesh::ExplicitSystem>(system_name);
                }
            }

            void describe(std::ostream &os) const override { os << "Sys"; }

            std::string system_name{"main"};
            std::string system_type{"linear_implicit"};

            std::vector<Var> vars;
            std::vector<BC> bcs;
        };

        FunctionSpace::FunctionSpace(const Communicator &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
            impl_->systems = std::make_shared<libMesh::EquationSystems>(mesh->raw_type());
        }

        FunctionSpace::~FunctionSpace() {}

        void FunctionSpace::read(Input &in) {
            Sys main_system;
            main_system.read(in);

            if (impl_->mesh->empty()) {
                in.get("mesh", *impl_->mesh);

                // Users did not care, so lets give them a cube
                if (impl_->mesh->empty()) {
                    impl_->mesh->unit_cube();
                }
            }

            if (!impl_->systems) {
                impl_->systems = std::make_shared<libMesh::EquationSystems>(impl_->mesh->raw_type());
            }

            if (!Options().parse(in)) {
                return;
            }
        }

        void FunctionSpace::describe(std::ostream &os) const {}

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }

        const Mesh &FunctionSpace::mesh() const { return *impl_->mesh; }

        Mesh &FunctionSpace::mesh() { return *impl_->mesh; }

    }  // namespace libmesh
}  // namespace utopia
