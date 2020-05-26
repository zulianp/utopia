#include "utopia_PoissonApp.hpp"

#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_ui.hpp"

#include "utopia_Flow.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_TensorInterpolate.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh_refinement.h"

namespace utopia {

    void PoissonApp::run(Input &in) {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<LibMeshFunctionSpace, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);

        bool no_solve = false;
        in.get("no-solve", no_solve);

        auto &V = space.space().subspace(0);
        auto u = trial(V);
        auto v = test(V);

        USparseMatrix A;
        UVector rhs, x, forcing_term;

        auto linear_form = inner(coeff(0.0), v) * dX;
        assemble(linear_form, rhs);

        std::string path;
        in.get("tensor-data-path", path);

        int dim = V.mesh().spatial_dimension();

        UIndexSet ghost_nodes;
        convert(V.dof_map().get_send_list(), ghost_nodes);

        UVector tensor_diffusivity;
        if (!path.empty() && read_tensor_data(path, local_size(rhs), dim * dim, tensor_diffusivity)) {
            build_tensor_ghosted(
                V.dof_map().n_local_dofs(), V.dof_map().n_dofs(), dim * dim, ghost_nodes, tensor_diffusivity);

            auto D = tensor_interpolate(V, tensor_diffusivity, dim, dim);
            auto bilinear_form = inner(D * grad(u), grad(v)) * dX;
            assemble(bilinear_form, A);
        } else {
            auto f = ctx_fun(lambda_fun([](const std::vector<double> &p) -> double {
                // if(p[1] < 0.5) {
                //     return 10.0;
                // } else {
                //     return 0.001;
                // }
                return 0.003325;
            }));

            auto bilinear_form = inner(f * grad(u), grad(v)) * dX;
            assemble(bilinear_form, A);
        }

        if (no_solve) return;

        x = local_zeros(local_size(rhs));

        forcing_term = local_zeros(local_size(rhs));
        forcing_function.eval(x, forcing_term);
        rhs += forcing_term;

        apply_boundary_conditions(V, A, rhs);

        bool use_mg = false;
        in.get("use-mg", use_mg);

        if (!use_mg) {
            Factorization<USparseMatrix, UVector> fact(MATSOLVERMUMPS, PCLU);
            fact.describe(std::cout);
            fact.solve(A, rhs, x);
        } else {
            auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
            auto smoother = std::make_shared<GaussSeidel<USparseMatrix, UVector>>();

            in.get("gs", *smoother);
            // auto smoother = std::make_shared<ProjectedGaussSeidel<USparseMatrix, UVector>>();
            // auto smoother = std::make_shared<ConjugateGradient<USparseMatrix, UVector, HOMEMADE>>();
            // linear_solver->verbose(true);
            SemiGeometricMultigrid mg(smoother, linear_solver);
            // mg.algebraic().rtol(1e-16);
            // mg.algebraic().atol(1e-16);
            // mg.algebraic().max_it(400);

            int n_levels = 3;
            int max_it = 80;
            bool verbose = true;
            bool solve_problem = true;
            bool write_op = false;

            in.get("vebose", verbose);
            in.get("n-levels", n_levels);
            in.get("max-it", max_it);
            in.get("solve-problem", solve_problem);
            in.get("write-op", write_op);

            in.get("multigrid", mg);
            mg.max_it(max_it);
            mg.verbose(verbose);

            mg.init(V.equation_system(), n_levels);
            mg.update(make_ref(A));

            if (solve_problem) {
                mg.apply(rhs, x);
            }

            if (write_op) {
                mg.algebraic().write("./");
            }
        }

        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}  // namespace utopia
