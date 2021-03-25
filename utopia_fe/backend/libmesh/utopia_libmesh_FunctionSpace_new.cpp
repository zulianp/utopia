#include "utopia_libmesh_FunctionSpace_new.hpp"

#include "utopia_Options.hpp"

#include "utopia_libmesh_ConvertTensor.hpp"

#include <cstdlib>

// All libmesh includes
#include "libmesh/boundary_info.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/namebased_io.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/parsed_function.h"
#include "libmesh/reference_counter.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/utility.h"

namespace utopia {
    namespace libmesh {

        static const int MAIN_SYSTEM_ID = 0;

        class FunctionSpaceWrapper {
        public:
            virtual ~FunctionSpaceWrapper() = default;

            std::shared_ptr<Mesh> mesh;
            std::shared_ptr<libMesh::EquationSystems> systems;
        };

        class FunctionSubspace::Impl {
        public:
            std::shared_ptr<FunctionSpaceWrapper> space;
            int system_id{MAIN_SYSTEM_ID};
            int subspace_id{0};
            int n_vars{1};
        };

        class FunctionSpace::Var : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("fe_family", fe_family, "The finite element family from enum libMesh::FEFamily.")
                         .add_option("order", order, "The finite element order from enum libMesh::Order.")
                         .add_option("name", name, "Name of the variable")
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
                    static const int N_SUFFIXES = 4;
                    static const char *suffixes[N_SUFFIXES] = {"_x", "_y", "_z", "_t"};

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
        };

        class FunctionSpace::BC : public Configurable, public Describable {
        public:
            void read(Input &in) override {
                if (!Options()
                         .add_option("name", name, "Name of boundary condition.")
                         .add_option("type", type, "Type of boundary condition.")
                         .add_option("var", var, "Variable number associated with boundary condition.")
                         .add_option("function_type", function_type, "Type of function: const|parsed.")
                         .add_option("value", value, "Value of function")
                         .add_option("side", side, "The sideset id")
                         .parse(in)) {
                    return;
                }
            }

            void add_to_system(libMesh::System &sys) {
                std::vector<unsigned int> vars(1);
                vars[0] = var;
                std::set<libMesh::boundary_id_type> bt;

                auto &dof_map = sys.get_dof_map();

                if (side == -1) {
                    side = sys.get_mesh().get_boundary_info().get_id_by_name(name);

                    if (side == libMesh::BoundaryInfo::invalid_id) {
                        assert(false);
                        Utopia::Abort();
                    }
                }

                bt.insert(side);

                if (type == "dirichlet") {
                    if (function_type == "const") {
                        dof_map.add_dirichlet_boundary(libMesh::DirichletBoundary(
                            bt, vars, libMesh::ConstFunction<libMesh::Real>(atof(value.c_str()))));
                    } else {
                        dof_map.add_dirichlet_boundary(
                            libMesh::DirichletBoundary(bt, vars, libMesh::ParsedFunction<libMesh::Real>(value)));
                    }
                }
            }

            void describe(std::ostream &os) const override {
                os << "name: " << name << ", ";
                os << "type: " << type << ", ";
                os << "function_type: " << function_type << ", ";
                os << "value: " << value << ", ";
                os << "var: " << var << ", ";
                os << "side: " << side << '\n';
            }

            std::string name{"bc"};
            std::string type{"dirichlet"};
            std::string function_type{"const"};
            std::string value{"0"};
            int var{0};
            int side{-1};
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
                libMesh::System *sys;
                if ("linear_implicit" == system_type) {
                    sys = &impl.systems->add_system<libMesh::LinearImplicitSystem>(system_name);
                } else if ("nonlinear_implicit" == system_type) {
                    sys = &impl.systems->add_system<libMesh::NonlinearImplicitSystem>(system_name);
                } else
                // if ("explicit" == system_type)
                {
                    sys = &impl.systems->add_system<libMesh::ExplicitSystem>(system_name);
                }

                for (auto &v : vars) {
                    v.add_to_system(*sys);
                }

                for (auto &b : bcs) {
                    b.add_to_system(*sys);
                }
            }

            void describe(std::ostream &os) const override {
                os << "system_name: " << system_name << '\n';
                os << "system_type: " << system_type << '\n';

                for (auto &v : vars) {
                    v.describe(os);
                }

                for (auto &b : bcs) {
                    b.describe(os);
                }
            }

            std::string system_name{"main"};
            std::string system_type{"linear_implicit"};

            std::vector<Var> vars;
            std::vector<BC> bcs;
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
            impl_->systems = std::make_shared<libMesh::EquationSystems>(mesh->raw_type());
        }

        bool FunctionSpace::empty() const { return !static_cast<bool>(impl_->systems); }

        FunctionSpace::~FunctionSpace() {}

        void FunctionSpace::read(Input &in) {
            Sys main_system;
            main_system.read(in);

            if (!Options().parse(in)) {
                return;
            }

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

            main_system.add_to_space(*impl_);
            impl_->systems->init();
        }

        bool FunctionSpace::read_with_state(Input &in, Vector &val) {
            Sys main_system;
            main_system.read(in);

            if (!Options().parse(in)) {
                return false;
            }

            Path path;
            in.get("mesh", [&](Input &in) { in.get("path", path); });

            int time_step = 1;
            in.get("time_step", time_step);

            auto ext = path.extension();

            if (ext != "e") {
                utopia::err() << "Could not read mesh at " << path.to_string()
                              << ", only exodus (.e) format supported\n";
                return false;
            }

            if (mesh().empty()) {
                mesh().init_distributed();
            }

            libMesh::ExodusII_IO io(mesh().raw_type());
            io.read(path.c_str());

            mesh().set_database(path);

            mesh().raw_type().prepare_for_use();

            if (!impl_->systems) {
                impl_->systems = std::make_shared<libMesh::EquationSystems>(impl_->mesh->raw_type());
            }

            main_system.add_to_space(*impl_);
            impl_->systems->init();

            auto &system = impl_->systems->get_system(system_id());

            for (auto &v : main_system.vars) {
                io.copy_nodal_solution(system, v.name, v.name, time_step);
            }

            convert(*system.solution, val);
            return !val.empty();
        }

        bool FunctionSpace::read(const Path &path,
                                 const std::vector<std::string> &var_names,
                                 Vector &val,
                                 const int time_step) {
            auto ext = path.extension();

            if (ext != "e") {
                utopia::err() << "Could not read mesh at " << path.to_string()
                              << ", only exodus (.e) format supported\n";
                return false;
            }

            if (mesh().empty()) {
                mesh().init_distributed();
            }

            libMesh::ExodusII_IO io(mesh().raw_type());
            io.read(path.c_str());

            mesh().set_database(path);

            mesh().raw_type().prepare_for_use();

            Sys main_system;
            for (auto &v_name : var_names) {
                Var var;
                var.name = v_name;
                main_system.vars.push_back(var);
            }

            if (!impl_->systems) {
                impl_->systems = std::make_shared<libMesh::EquationSystems>(impl_->mesh->raw_type());
            }

            main_system.add_to_space(*impl_);
            impl_->systems->init();

            auto &system = impl_->systems->get_system(system_id());

            for (auto &v_name : var_names) {
                io.copy_nodal_solution(system, v_name, v_name, 1);
            }

            convert(*system.solution, val);
            return !val.empty();
        }

        void FunctionSpace::write(const Path &path, const Vector &x) {
            assert(impl_->systems);

            auto &sys = impl_->systems->get_system(system_id());
            convert(x, *sys.solution);

            // if (path.extension() == "e") {
            //     libMesh::ExodusII_IO(mesh().raw_type()).write_equation_systems(path.to_string(), *impl_->systems);
            // } else {
            libMesh::NameBasedIO(impl_->mesh->raw_type()).write_equation_systems(path.to_string(), *impl_->systems);
            // }
        }

        void FunctionSpace::describe(std::ostream &os) const {
            if (!mesh().empty()) {
                mesh().describe(os);
            }
        }

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }

        const Mesh &FunctionSpace::mesh() const { return *impl_->mesh; }

        Mesh &FunctionSpace::mesh() { return *impl_->mesh; }

        libMesh::EquationSystems &FunctionSpace::raw_type() {
            assert(impl_->systems);
            return *impl_->systems;
        }

        const libMesh::EquationSystems &FunctionSpace::raw_type() const {
            assert(impl_->systems);
            return *impl_->systems;
        }

        const libMesh::DofMap &FunctionSpace::raw_type_dof_map() const {
            return impl_->systems->get_system(system_id()).get_dof_map();
        }

        void FunctionSpace::create_matrix(Matrix &mat) const {
            UTOPIA_TRACE_REGION_BEGIN("libmesh::FunctionSpace::create_matrix");

            using IndexSet = Traits<Vector>::IndexSet;

            auto &sys = impl_->systems->get_system(system_id());
            auto &dof_map = sys.get_dof_map();

            auto vl = layout(comm(), dof_map.n_local_dofs(), dof_map.n_local_dofs());
            auto ml = square_matrix_layout(vl);

            IndexSet d_nnz, o_nnz;
            convert(dof_map.get_n_nz(), d_nnz);
            convert(dof_map.get_n_oz(), o_nnz);

            mat.sparse(ml, d_nnz, o_nnz);

            UTOPIA_TRACE_REGION_END("libmesh::FunctionSpace::create_matrix");
        }

        void FunctionSpace::create_vector(Vector &vec) const {
            UTOPIA_TRACE_REGION_BEGIN("libmesh::FunctionSpace::create_vector");

            using IndexSet = Traits<Vector>::IndexSet;
            auto &sys = impl_->systems->get_system(system_id());
            auto &dof_map = sys.get_dof_map();

            IndexSet index;
            convert(dof_map.get_send_list(), index);

            auto vl = layout(comm(), dof_map.n_local_dofs(), dof_map.n_dofs());
            vec.ghosted(vl, index);

            UTOPIA_TRACE_REGION_END("libmesh::FunctionSpace::create_vector");
        }

        void FunctionSpace::apply_constraints(Matrix &mat, Vector &vec) const {
            UTOPIA_TRACE_REGION_BEGIN("libmesh::FunctionSpace::apply_constraints");

            assert(!utopia::empty(mat));
            assert(!utopia::empty(vec));

            auto &dof_map = impl_->systems->get_system(system_id()).get_dof_map();

            const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

            Size ls = local_size(mat);
            Size s = size(mat);

            Traits<Matrix>::IndexSet index;

            Range rr = range(vec);

            if (has_constaints) {
                for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                    if (dof_map.is_constrained_dof(i)) {
                        index.push_back(i);
                    }
                }
            }

            set_zero_rows(mat, index, 1.);

            Write<Vector> w_v(vec);

            if (has_constaints) {
                libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

                Range r = range(vec);
                for (SizeType i = r.begin(); i < r.end(); ++i) {
                    if (dof_map.is_constrained_dof(i)) {
                        auto valpos = rhs_values.find(i);
                        vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                    }
                }
            }

            UTOPIA_TRACE_REGION_END("libmesh::FunctionSpace::apply_constraints");
        }

        FunctionSpace::SizeType FunctionSpace::system_id() const { return MAIN_SYSTEM_ID; }

        std::shared_ptr<FunctionSpaceWrapper> FunctionSpace::wrapper() { return impl_; }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const {
            assert(impl_->systems);
            return impl_->systems->get_system(system_id()).n_dofs();
        }

        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {
            assert(impl_->systems);
            return impl_->systems->get_system(system_id()).n_local_dofs();
        }

        FunctionSpace::SizeType FunctionSpace::n_subspaces() const {
            assert(impl_->systems);
            return impl_->systems->get_system(system_id()).n_vars();
        }

        FunctionSpace::SizeType FunctionSpace::n_systems() const {
            assert(impl_->systems);
            return impl_->systems->n_systems();
        }

        FunctionSubspace FunctionSpace::subspace(const SizeType i, const SizeType n_vars) {
            assert(i < n_subspaces());
            assert(i + n_vars <= n_subspaces());
            FunctionSubspace ret;
            ret.impl_->system_id = system_id();
            ret.impl_->subspace_id = i;
            ret.impl_->n_vars = n_vars;
            ret.impl_->space = wrapper();
            return ret;
        }

        FunctionSubspace FunctionSpace::auxiliary_space(const SizeType i) {
            FunctionSubspace ret;
            assert(i < n_systems());

            ret.impl_->system_id = i;
            ret.impl_->subspace_id = 0;
            ret.impl_->space = wrapper();
            return ret;
        }

        FunctionSubspace::FunctionSubspace() : impl_(std::make_shared<Impl>()) {}

        FunctionSubspace::~FunctionSubspace() {}

        FunctionSubspace::SizeType FunctionSubspace::subspace_id() const { return impl_->subspace_id; }

        FunctionSubspace::SizeType FunctionSubspace::n_dofs() const {
            return impl_->space->systems->get_system(impl_->system_id).n_dofs();
        }

        FunctionSubspace::SizeType FunctionSubspace::n_local_dofs() const {
            return impl_->space->systems->get_system(impl_->system_id).n_local_dofs();
        }

        FunctionSubspace::SizeType FunctionSubspace::n_subspaces() const { return impl_->n_vars; }

        // access main function space subspaces (main system)
        FunctionSubspace FunctionSubspace::subspace(const SizeType i, const SizeType n_vars) {
            assert(i < n_subspaces());
            assert(i + n_vars <= n_subspaces());
            FunctionSubspace ret;
            ret.impl_->subspace_id = i;
            ret.impl_->n_vars = n_vars;
            ret.impl_->space = this->impl_->space;
            return ret;
        }

    }  // namespace libmesh
}  // namespace utopia
