#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_Options.hpp"

// All libmesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/reference_counter.h"

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
                         .parse(in)) {
                    return;
                }
            }

            void describe(std::ostream &os) const override { os << "Var"; }

            std::string fe_family{"LAGRANGE"};
            std::string order{"FIRST"};
            std::string name{"var"};
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

            void describe(std::ostream &os) const override { os << "BC"; }
        };

        class FunctionSpace::Sys : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("system_name", system_name, "Name of system.")
                         .add_option("system_type", system_type, "Type of system associated with function space.")
                         .parse(in)) {
                    return;
                }
            }

            void describe(std::ostream &os) const override { os << "Sys"; }

            std::string system_name{"sys"};
            std::string system_type{"linear_implicit"};

            std::vector<Var> vars;
            std::vector<BC> bcs;
        };

        FunctionSpace::FunctionSpace(const Communicator &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
        }

        FunctionSpace::~FunctionSpace() {}

        void FunctionSpace::read(Input &in) {
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
