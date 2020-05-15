#ifndef UTOPIA_LIBMESH_NON_LINEAR_FE_FUNCTION_HPP
#define UTOPIA_LIBMESH_NON_LINEAR_FE_FUNCTION_HPP

#include "utopia_NonLinearFEFunction.hpp"
#include "utopia_fe_base.hpp"

#include "utopia.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/parallel_mesh.h"
#include "utopia_libmesh.hpp"

#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_Projection.hpp"
#include "utopia_libmesh_Assembler.hpp"
#include "utopia_libmesh_FEBackend.hpp"

namespace utopia {
    // libmesh
    template <class FunctionSpaceT, class Expr, class GlobalMatrix, class GlobalVector>
    void element_assemble_expression_v(const libMesh::MeshBase::const_element_iterator &it,
                                       const Expr &expr,
                                       GlobalMatrix &mat,
                                       GlobalVector &vec,
                                       const bool apply_constraints = true) {
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        static const int Backend = TraitsT::Backend;

        const auto &space = find_space<FunctionSpaceT>(expr);
        const auto &dof_map = space.dof_map();

        AssemblyContext<Backend> ctx;
        ctx.set_current_element((*it)->id());
        ctx.set_has_assembled(false);

        ElementMatrix el_mat;
        ElementVector el_vec;

        ctx.init_bilinear(expr);

        FormEvaluator<Backend> eval;
        eval.eval(expr, el_mat, el_vec, ctx);

        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*it, dof_indices);

        if (ctx.has_assembled()) {
            // if(apply_constraints) {
            //     dof_map.heterogenously_constrain_element_matrix_and_vector(el_mat.implementation(),
            //     el_vec.implementation(), dof_indices);
            // }

            add_matrix(el_mat, dof_indices, dof_indices, mat);
            add_vector(el_vec, dof_indices, vec);
        }
    }

    // libmesh
    template <class FunctionSpaceT, class Expr, class GlobalMatrix>
    void element_assemble_expression_v(const libMesh::MeshBase::const_element_iterator &it,
                                       const Expr &expr,
                                       Tensor<GlobalMatrix, 2> &t_mat) {
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        static const int Backend = TraitsT::Backend;

        auto &mat = t_mat.derived();
        const auto &space = find_space<FunctionSpaceT>(expr);
        const auto &dof_map = space.dof_map();

        AssemblyContext<Backend> ctx;
        ctx.set_current_element((*it)->id());
        ctx.set_has_assembled(false);

        ElementMatrix el_mat;
        ctx.init_bilinear(expr);

        FormEvaluator<Backend> eval;
        eval.eval(expr, el_mat, ctx, true);

        if (ctx.has_assembled()) {
            std::vector<libMesh::dof_id_type> dof_indices;
            dof_map.dof_indices(*it, dof_indices);

            add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
        }
    }

    template <class FunctionSpaceT, class Expr, typename T>
    void element_assemble_expression_v(const libMesh::MeshBase::const_element_iterator &it,
                                       const Expr &expr,
                                       libMesh::SparseMatrix<T> &mat) {
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        static const int Backend = TraitsT::Backend;

        const auto &space = find_space<FunctionSpaceT>(expr);
        const auto &dof_map = space.dof_map();

        AssemblyContext<Backend> ctx;
        ctx.set_current_element((*it)->id());
        ctx.set_has_assembled(false);

        ElementMatrix el_mat;
        ctx.init_bilinear(expr);

        FormEvaluator<Backend> eval;
        eval.eval(expr, el_mat, ctx, true);

        if (ctx.has_assembled()) {
            std::vector<libMesh::dof_id_type> dof_indices;
            dof_map.dof_indices(*it, dof_indices);

            add_matrix(el_mat.implementation(), dof_indices, dof_indices, mat);
        }
    }

    template <class FunctionSpaceT, class Expr, class GlobalVector>
    void element_assemble_expression_v(const libMesh::MeshBase::const_element_iterator &it,
                                       const Expr &expr,
                                       Tensor<GlobalVector, 1> &t_vec) {
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;

        static const int Backend = TraitsT::Backend;

        const auto &space = find_space<FunctionSpaceT>(expr);
        const auto &dof_map = space.dof_map();
        auto &vec = t_vec.derived();

        AssemblyContext<Backend> ctx;
        ctx.set_current_element((*it)->id());

        ElementVector el_vec;

        ctx.init_linear(expr);

        FormEvaluator<Backend> eval;
        eval.eval(expr, el_vec, ctx, true);

        std::vector<libMesh::dof_id_type> dof_indices;
        dof_map.dof_indices(*it, dof_indices);

        add_vector(el_vec.implementation(), dof_indices, vec);
    }

    // homemade
    template <class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
    void element_assemble_expression_v(const int element_index,
                                       const Equality<Left, Right> &equation,
                                       Matrix &mat,
                                       Vector &rhs,
                                       const bool apply_constraints = true) {
        std::cout << element_index << std::endl;
    }

    template <class FunctionSpaceT, class Left, class Right, class Matrix, class Vector>
    void assemble_expression_v(const Equality<Left, Right> &equation,
                               Matrix &mat,
                               Vector &vec,
                               const bool apply_constraints = true) {
        auto &space = find_space<FunctionSpaceT>(equation);
        space.initialize();

        auto &m = space.mesh();
        auto &dof_map = space.dof_map();

        auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
                                  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

        mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
        vec = local_zeros(dof_map.n_local_dofs());

        Write<Matrix> w_m(mat);
        Write<Vector> w_v(vec);

        for (auto it = elements_begin(m); it != elements_end(m); ++it) {
            element_assemble_expression_v<FunctionSpaceT>(it, equation, mat, vec, apply_constraints);
        }
    }

    template <class Expr>
    bool assemble(const Expr &expr, double &val) {
        return LibMeshAssembler().assemble(expr, val);
    }

    template <class Expr>
    bool assemble(Expr &expr, double &val) {
        return LibMeshAssembler().assemble(expr, val);
    }

    template <class Expr>
    bool assemble(const Expr &expr, USparseMatrix &mat, const bool first = true) {
        return LibMeshAssembler().assemble(expr, mat);
    }

    template <class Expr>
    bool assemble(Expr &expr, USparseMatrix &mat, const bool first = true) {
        return LibMeshAssembler().assemble(expr, mat);
    }

    template <class Expr, typename T>
    bool assemble(const Expr &expr, libMesh::SparseMatrix<T> &mat, const bool first = true) {
        typedef typename FindFunctionSpace<Expr>::Type FunctionSpaceT;
        auto &space = find_space<FunctionSpaceT>(expr);

        if (!first) {
            mat.zero();
        }

        {
            auto &m = space.mesh();
            for (auto it = elements_begin(m); it != elements_end(m); ++it) {
                element_assemble_expression_v<FunctionSpaceT>(it, expr, mat);
            }
        }

        return true;
    }

    template <class Expr>
    bool assemble(const Expr &expr, UVector &vec, const bool first = true) {
        return LibMeshAssembler().assemble(expr, vec);
    }

    template <class Expr>
    bool assemble(Expr &expr, UVector &vec, const bool first = true) {
        return LibMeshAssembler().assemble(expr, vec);
    }

    template <class... Eqs>
    bool assemble(const Equations<Eqs...> &eqs, USparseMatrix &mat, UVector &vec) {
        return LibMeshAssembler().assemble(eqs, mat, vec);
    }

    template <class Left, class Right>
    bool assemble(const Equality<Left, Right> &equation, USparseMatrix &mat, UVector &vec) {
        return assemble(equations(equation), mat, vec);
    }

    template <class FunctionSpaceT>
    bool assemble(EquationIntegrator<FunctionSpaceT> &equation, USparseMatrix &mat, UVector &vec) {
        return LibMeshAssembler().assemble(equation, mat, vec);
    }

    template <class FunctionSpaceT>
    bool assemble(const EquationIntegrator<FunctionSpaceT> &equation, USparseMatrix &mat, UVector &vec) {
        return LibMeshAssembler().assemble(equation, mat, vec);
    }

    template <class Expr>
    void init_constraints(const Expr &expr) {
        FEBackend<LIBMESH_TAG>::init_constraints(expr);
    }

    template <class From, class Tensor>
    void apply(const Projection<From, Tensor> &expr) {
        typedef typename FindFunctionSpace<From>::Type FunctionSpaceT;
        auto &space = find_space<FunctionSpaceT>(expr.from());

        // init vector
        // assemble weighted coeffs
        // remove mass contrib
    }

    template <class Matrix, class Vector, class Eqs>
    class NonLinearFEFunction : public Function<Matrix, Vector> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        NonLinearFEFunction(const Eqs &eqs, const bool compute_linear_residual = false)
            : eqs_(eqs), first_(true), compute_linear_residual_(compute_linear_residual) {}

        virtual ~NonLinearFEFunction() {}

        virtual bool value(const Vector &x, Scalar &value) const override {
            value = dot(x, buff_mat * x) - dot(x, buff_vec);
            return true;
        }

        virtual bool gradient(const Vector &x, Vector &result) const override {
            if (compute_linear_residual_) {
                result = buff_mat * x - buff_vec;
            } else {
                result = buff_vec;
            }

            return true;
        }

        virtual bool hessian(const Vector &, Matrix &H) const override {
            H = buff_mat;
            return true;
        }

        virtual bool update(const Vector &x) override {
            typedef decltype(eqs_.template get<0>()) Eq1;
            typedef typename FindFunctionSpace<Eq1>::Type FunctionSpaceT;
            auto &space = find_space<FunctionSpaceT>(eqs_);

            // FIXME
            auto &x_mutable = const_cast<Vector &>(x);
            synchronize(x_mutable);  //.implementation().update_ghosts();

            LibMeshAssembler().assemble(eqs_, buff_mat, buff_vec);

            // if(first_) {
            // 	auto &dof_map = space.dof_map();
            // 	auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
            // 		*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

            // 	buff_mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
            // 	buff_vec = local_zeros(dof_map.n_local_dofs());
            // } else {
            // 	buff_mat *= 0.;
            // 	buff_vec *= 0.;
            // }

            // {
            // 	Write<USparseMatrix> w_m(buff_mat);
            // 	Write<UVector>  w_v(buff_vec);

            // 	auto &m = space.mesh();

            // 	for(auto it = elements_begin(m); it != elements_end(m); ++it) {
            // 		element_assemble_expression_v<FunctionSpaceT>(it, eqs_, buff_mat, buff_vec, false);
            // 	}
            // }

            if (first_ && !compute_linear_residual_) {
                buff_vec *= -1.;
                apply_boundary_conditions(space.dof_map(), buff_mat, buff_vec);
                buff_vec *= -1.;
            } else {
                apply_boundary_conditions(space.dof_map(), buff_mat, buff_vec);
            }

            if (!first_ && !compute_linear_residual_) {
                apply_zero_boundary_conditions(space.dof_map(), buff_vec);
            }

            first_ = false;
            return true;
        }

        void reset() { first_ = true; }

        Eqs eqs_;
        bool first_;

        Matrix buff_mat;
        Vector buff_vec;

        bool compute_linear_residual_;
    };

    template <class... Eqs, class... Constr>
    bool nl_solve(const Equations<Eqs...> &eqs, const FEConstraints<Constr...> &constr, UVector &sol) {
        typedef typename GetFirst<Eqs...>::Type Eq1Type;
        typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;

        FEBackend<LIBMESH_TAG>::init_constraints(constr);

        auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());
        space.initialize();

        UIndexSet ghost_nodes;
        convert(space.dof_map().get_send_list(), ghost_nodes);
        sol = ghosted(space.dof_map().n_local_dofs(), space.dof_map().n_dofs(), ghost_nodes);

        NonLinearFEFunction<USparseMatrix, UVector, Equations<Eqs...>> nl_fun(eqs);
        Newton<USparseMatrix, UVector> solver(std::make_shared<Factorization<USparseMatrix, UVector>>());
        // solver.set_line_search_strategy(std::make_shared<Backtracking<USparseMatrix, UVector>>());
        solver.verbose(true);
        return solver.solve(nl_fun, sol);
    }

    template <class... Eqs, class... Constr>
    bool nl_implicit_euler(const Equations<Eqs...> &eqs,
                           const FEConstraints<Constr...> &constr,
                           UVector &old_sol,
                           UVector &sol,
                           const double dt,
                           std::size_t n_ts = 10) {
        typedef typename GetFirst<Eqs...>::Type Eq1Type;
        typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;

        FEBackend<LIBMESH_TAG>::init_constraints(constr);

        auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());
        space.initialize();
        auto &dof_map = space.dof_map();

        UIndexSet ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);

        sol = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
        old_sol = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);

        NonLinearFEFunction<USparseMatrix, UVector, Equations<Eqs...>> nl_fun(eqs, true);
        Newton<USparseMatrix, UVector> solver(std::make_shared<Factorization<USparseMatrix, UVector>>());
        solver.verbose(true);

        libMesh::ExodusII_IO io(space.mesh());

        for (std::size_t ts = 0; ts < n_ts; ++ts) {
            nl_fun.reset();
            if (!solver.solve(nl_fun, sol)) return false;
            old_sol = sol;

            convert(sol, *space.equation_system().solution);
            space.equation_system().solution->close();
            io.write_timestep(space.equation_system().name() + ".e", space.equation_systems(), ts + 1, ts * dt);
        }

        return true;
    }

    template <class... Eqs, class... Constr>
    bool solve(const Equations<Eqs...> &eqs,
               const FEConstraints<Constr...> &constr,
               UVector &sol,
               const bool first = true) {
        // FIXME this stuff only works for the libmesh backend
        typedef typename GetFirst<Eqs...>::Type Eq1Type;
        typedef typename FindFunctionSpace<Eq1Type>::Type FunctionSpaceT;
        auto &space = find_space<FunctionSpaceT>(eqs.template get<0>());

        init_constraints(constr);

        space.initialize();
        auto &m = space.mesh();
        auto &dof_map = space.dof_map();
        auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
                                  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

        USparseMatrix mat;
        UVector vec;

        mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
        vec = local_zeros(dof_map.n_local_dofs());

        if (empty(sol)) {
            sol = local_zeros(local_size(vec));
            std::cout << "[Warning] empty solution vector initializing to zero" << std::endl;
        }

        {
            Write<USparseMatrix> w_m(mat);
            Write<UVector> w_v(vec);
            Read<UVector> r_s(sol);

            for (auto it = elements_begin(m); it != elements_end(m); ++it) {
                element_assemble_expression_v<FunctionSpaceT>(it, eqs, mat, vec, false);
            }
        }

        apply_boundary_conditions(dof_map, mat, vec);

        if (!first) {
            apply_zero_boundary_conditions(dof_map, vec);
        }

        UVector inc = local_zeros(local_size(sol));
        Factorization<USparseMatrix, UVector> solver;
        if (!solver.solve(mat, vec, inc)) {
            return false;
        }

        sol += inc;

        double norm_r = norm2(inc);
        std::cout << "norm(inc): " << norm_r << std::endl;
        return true;
    }

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_NON_LINEAR_FE_FUNCTION_HPP
