#include "utopia_FlowWIthFracturesApp.hpp"

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
#include "utopia_FlowWithFractures.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"

namespace utopia {

    void FlowWithFracturesApp::run(Input &in)
    {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        auto &V = space.space().subspace(0);

        FlowWithFractures<LibMeshFunctionSpace, USparseMatrix, UVector> model(V);
    
        in.get("model", model);

        UVector x = local_zeros(V.dof_map().n_local_dofs()), rhs;
        USparseMatrix A;
        model.assemble_hessian_and_gradient(x, A, rhs);
        apply_boundary_conditions(V, A, rhs);
      
        Factorization<USparseMatrix, UVector> fact(MATSOLVERMUMPS,PCLU);
        fact.describe(std::cout);
        fact.solve(A, rhs, x);

        rename("a", A);
        write("A.m", A);

        rename("b", rhs);
        write("B.m", rhs);

        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}
