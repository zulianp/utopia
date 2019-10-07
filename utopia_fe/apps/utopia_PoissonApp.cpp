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
        


        auto &V = space.space().subspace(0);
       

        auto u = trial(V);
        auto v = test(V);

        auto linear_form = inner(coeff(0.0), v) * dX;
        auto bilinear_form = inner(grad(u), grad(v)) * dX;

        USparseMatrix A;
        UVector rhs, x, neumann_bc;
        

        //rhs.set(1.0);
        assemble(bilinear_form == linear_form, A, rhs);

        x = local_zeros(local_size(rhs));

        neumann_bc = local_zeros(local_size(rhs));
        
        forcing_function.eval(x,neumann_bc);

        rhs+=neumann_bc;

        // utopia::disp(size(A).get(0));

        // utopia::disp(size(A).get(1));

        // utopia::write("A_before.m", A);

        apply_boundary_conditions(V, A, rhs);


        // utopia::rename("a", A);

        // utopia::write("A.m", A);

        // utopia::rename("b", rhs);

        // utopia::write("rhs.m", rhs);

      

        Factorization<USparseMatrix, UVector> fact(MATSOLVERMUMPS,PCLU);
        fact.describe(std::cout);
        fact.solve(A, rhs, x);


        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}
