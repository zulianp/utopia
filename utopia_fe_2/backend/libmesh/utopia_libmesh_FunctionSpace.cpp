#include "utopia_libmesh_FunctionSpace.hpp"

// utopia
#include "utopia_Readable.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Writable.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"
#include "utopia_petsc_Matrix.hpp"

// libmesh
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/string_to_enum.h"

namespace utopia {

    template <class Vector>
    inline void convert(const Vector &utopia_vec, libMesh::NumericVector<libMesh::Number> &lm_vec) {
        Read<Vector> w_s(utopia_vec);
        Range r = range(utopia_vec);
        for (long i = r.begin(); i < r.end(); ++i) {
            lm_vec.set(i, utopia_vec.get(i));
        }
    }

    void FunctionSpace<LMMesh>::create_vector(Vector &x) const {
        auto &sys = this->system();
        auto &dof_map = sys.get_dof_map();

        Communicator comm(mesh_->comm().get());
        Traits::IndexArray ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);
        x.ghosted(layout(comm, dof_map.n_local_dofs(), dof_map.n_dofs()), ghost_nodes);
    }

    void FunctionSpace<LMMesh>::create_matrix(Matrix &A) const {
        auto &sys = this->system();
        auto &dof_map = sys.get_dof_map();

        Communicator comm(mesh_->comm().get());

        auto &lm_d_nnz = dof_map.get_n_nz();
        auto &lm_o_nnz = dof_map.get_n_oz();

        Traits::IndexArray d_nnz, o_nnz;

        convert(lm_d_nnz, d_nnz);
        convert(lm_o_nnz, o_nnz);

        auto l = layout(comm, dof_map.n_local_dofs(), dof_map.n_dofs());

        A.sparse(square_matrix_layout(l), d_nnz, o_nnz);
    }

    void FunctionSpace<LMMesh>::add_matrix(const Elem &e, const ElementMatrix &el_mat, Matrix &mat) const {
        auto &sys = this->system();
        auto &dof_map = sys.get_dof_map();

        std::vector<libMesh::dof_id_type> dofs;
        dof_map.dof_indices(&e, dofs);

        el_mat.read([&](const SizeType &i, const SizeType &j, const Scalar &val) { mat.c_add(dofs[i], dofs[j], val); });
    }

    void FunctionSpace<LMMesh>::add_vector(const Elem &e, const ElementVector &el_vec, Vector &vec) const {
        auto &sys = this->system();
        auto &dof_map = sys.get_dof_map();

        std::vector<libMesh::dof_id_type> dofs;
        dof_map.dof_indices(&e, dofs);

        const SizeType n = el_vec.size();

        assert(n == SizeType(dofs.size()));
        for (SizeType i = 0; i < n; ++i) {
            vec.c_add(dofs[i], el_vec.get(i));
        }
    }

    void FunctionSpace<LMMesh>::apply_constraints(Vector &vec) const {
        auto &dof_map = system().get_dof_map();

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<Vector> w_v(vec);

        if (has_constaints) {
            const libMesh::DofConstraintValueMap &rhs_values =
                const_cast<libMesh::DofMap &>(dof_map).get_primal_constraint_values();

            auto r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }
    }

    void FunctionSpace<LMMesh>::apply_constraints(Matrix &mat, Vector &vec) const {
        assert(!empty(mat));
        assert(!empty(vec));

        auto &dof_map = system().get_dof_map();

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        auto ls = local_size(mat);
        auto s = size(mat);

        Traits::IndexArray index;

        auto rr = range(vec);

        if (has_constaints) {
            for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    index.push_back(i);
                }
            }
        }

        mat.set_zero_rows(index, 1.);

        Write<Vector> w_v(vec);

        if (has_constaints) {
            const libMesh::DofConstraintValueMap &rhs_values =
                const_cast<libMesh::DofMap &>(dof_map).get_primal_constraint_values();

            auto r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }
    }

    void FunctionSpace<LMMesh>::apply_zero_constraints(Vector &vec) const {
        auto &dof_map = system().get_dof_map();

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<Vector> w_v(vec);

        if (has_constaints) {
            auto r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    vec.set(i, 0.0);
                }
            }
        }
    }

    class FunctionSpace<LMMesh>::Impl {
    public:
        int system_num;
        std::vector<int> vars;
    };

    libMesh::System &FunctionSpace<LMMesh>::system() { return equation_systems_->get_system(impl_->system_num); }

    const libMesh::System &FunctionSpace<LMMesh>::system() const {
        return equation_systems_->get_system(impl_->system_num);
    }

    FunctionSpace<LMMesh>::FunctionSpace(const Communicator &comm)
        : mesh_(std::make_shared<LMMesh>(comm)), impl_(utopia::make_unique<Impl>()) {}

    FunctionSpace<LMMesh>::~FunctionSpace() {}

    void FunctionSpace<LMMesh>::describe(std::ostream &os) const { mesh_->describe(os); }

    bool FunctionSpace<LMMesh>::write(const Path &path, const Vector &x) const {
        auto &sys = this->system();

        utopia::convert(x, *sys.solution);

        sys.solution->close();

        if (mesh_->comm().size() == 1) {
            libMesh::ExodusII_IO(mesh_->raw_type()).write_equation_systems(path.to_string(), *equation_systems_);
        } else {
            libMesh::Nemesis_IO(mesh_->raw_type()).write_equation_systems(path.to_string(), *equation_systems_);
        }

        return true;
    }

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
