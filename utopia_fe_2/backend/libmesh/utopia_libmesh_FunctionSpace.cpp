#include "utopia_libmesh_FunctionSpace.hpp"

// utopia
#include "utopia_SymbolicFunction.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"

// libmesh
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/string_to_enum.h"

namespace utopia {
    class FunctionSpace<LMMesh>::Impl {
    public:
        int system_num;
        std::vector<int> vars;
    };

    FunctionSpace<LMMesh>::FunctionSpace(const Communicator &comm)
        : mesh_(std::make_shared<LMMesh>(comm)), impl_(utopia::make_unique<Impl>()) {}

    FunctionSpace<LMMesh>::~FunctionSpace() {}

    void FunctionSpace<LMMesh>::describe(std::ostream &os) const { mesh_->describe(os); }

    void FunctionSpace<LMMesh>::read(Input &in) {
        mesh_->read(in);

        int n_vars = 0;

        std::vector<std::string> var_names;
        std::vector<libMesh::Order> var_orders;
        std::vector<libMesh::FEFamily> fe_families;

        std::string system_name = "main";
        libMesh::FEFamily fe_family = libMesh::LAGRANGE;
        libMesh::Order order = libMesh::FIRST;

        in.get("variables", [&](Input &sub_in_mid) {
            sub_in_mid.get_all([&](Input &sub_in) {
                std::string var_name = "u_" + std::to_string(n_vars);
                std::string var_order = "FIRST";
                std::string var_fe_family = "LAGRANGE";

                sub_in.get("name", var_name);
                sub_in.get("order", var_order);
                sub_in.get("fe_family", var_fe_family);

                var_names.push_back(var_name);
                var_orders.push_back(libMesh::Utility::string_to_enum<libMesh::Order>(var_order));
                fe_families.push_back(libMesh::Utility::string_to_enum<libMesh::FEFamily>(var_fe_family));

                ++n_vars;
            });
        });

        if (n_vars == 0) {
            var_names.push_back("u");
            var_orders.push_back(order);
            fe_families.push_back(fe_family);
            n_vars = 1;
        }

        in.get("system_name", system_name);

        if (!equation_systems_) {
            equation_systems_ = std::make_shared<libMesh::EquationSystems>(mesh_->raw_type());
        }

        // FIXME (parametrize system)
        auto &sys = equation_systems_->add_system<libMesh::LinearImplicitSystem>(system_name);
        impl_->system_num = sys.number();
        impl_->vars.resize(n_vars);

        for (int i = 0; i < n_vars; ++i) {
            const int var_num = sys.add_variable(var_names[i], var_orders[i], fe_families[i]);
            impl_->vars[i] = var_num;
        }

        in.get("boundary_conditions", [&](Input &in) {
            in.get_all([&](Input &in) {
                int side_set = 0;

                in.get("side", side_set);

                double value = 0;
                in.get("value", value);

                int var_num = 0;

                in.get("var", var_num);

                std::string type = "constant";  // expr

                in.get("type", type);

                if (type == "expr") {
#ifdef WITH_TINY_EXPR
                    std::string expr = "0";
                    in.get("value", expr);
                    auto f = symbolic(expr);

                    std::function<libMesh::Real(const libMesh::Point &)> fun =
                        [f](const libMesh::Point &xyz) -> libMesh::Real {
                        // FIXME
                        auto f_copy = f;
                        return f_copy.eval(xyz(0), xyz(1), xyz(2));
                    };

                    LibMeshLambdaFunction<libMesh::Real> lambda(fun);

                    std::vector<unsigned int> vars(1);
                    vars[0] = var_num;
                    std::set<libMesh::boundary_id_type> bt;
                    bt.insert(side_set);

                    libMesh::DirichletBoundary d_bc(bt, vars, lambda);
                    sys.get_dof_map().add_dirichlet_boundary(d_bc);

#else
                    assert(false);
                    std::cerr << "needs tiny expr library" << std::endl;
#endif

                } else {
                    std::vector<unsigned int> vars(1);
                    vars[0] = var_num;
                    std::set<libMesh::boundary_id_type> bt;
                    bt.insert(side_set);
                    libMesh::DirichletBoundary d_bc(bt, vars, libMesh::ConstFunction<libMesh::Real>(value));
                    sys.get_dof_map().add_dirichlet_boundary(d_bc);
                }
            });
        });

        sys.init();
    }

}  // namespace utopia
