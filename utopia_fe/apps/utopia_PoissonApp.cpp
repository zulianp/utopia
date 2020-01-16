#include "utopia_PoissonApp.hpp"

#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_Flow.hpp"
#include "utopia_TensorInterpolate.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"

namespace utopia {

    void PoissonApp::run(Input &in)
    {
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
        if(!path.empty() && read_tensor_data(path, local_size(rhs), dim*dim, tensor_diffusivity))
        {
            build_tensor_ghosted(
                V.dof_map().n_local_dofs(),
                V.dof_map().n_dofs(),
                dim*dim,
                ghost_nodes,
                tensor_diffusivity
            );

            auto D = tensor_interpolate(V, tensor_diffusivity, dim, dim);
            auto bilinear_form = inner(D * grad(u), grad(v)) * dX;
            assemble(bilinear_form, A);
        } else {
            auto bilinear_form = inner(grad(u), grad(v)) * dX;
            assemble(bilinear_form, A);
        }

        if(no_solve) return;

        x = local_zeros(local_size(rhs));

        forcing_term = local_zeros(local_size(rhs));
        forcing_function.eval(x, forcing_term);
        rhs += forcing_term;

        apply_boundary_conditions(V, A, rhs);

        Factorization<USparseMatrix, UVector> fact(MATSOLVERMUMPS,PCLU);
        fact.describe(std::cout);
        fact.solve(A, rhs, x);


        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}
