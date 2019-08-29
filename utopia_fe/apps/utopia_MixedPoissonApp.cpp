#include "utopia_MixedPoissonApp.hpp"

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

    void MixedPoissonApp::run(Input &in)
    {
        using FS   = utopia::LibMeshFunctionSpace;
        using PFS = utopia::ProductFunctionSpace<FS>;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<FS> space(make_ref(mesh));

        in.get("mesh", mesh);
        in.get("space", space);

   

        auto dim = mesh.mesh().mesh_dimension();
        auto &M  = space.space();


        UIForcingFunction<PFS, UVector> forcing_function(M);
        in.get("forcing-functions", forcing_function);
        
        if(dim + 1 != M.n_subspaces()) {
            std::cerr << "[Error] wrong number of subspaces, expected " <<   (dim + 1) << std::endl;
            std::cerr << "V[0] for scalar, V[1, ..., dim + 1] for gradient" << std::endl;
            return;
        }

        auto &V = M.subspace(0);
        auto W  = M.subspace(1, dim + 1);

        auto p = trial(V);
        auto q = test(V);

        auto v = trial(W);
        auto w = test(W);
        

        // does not work for some reason
        // auto bilinear_form = (inner(v, w)      * dX) -
        //                      (inner(p, div(w)) * dX) +
        //                      (inner(div(v), q) * dX);

        auto bilinear_form = (inner(v, w)      * dX)  +
                             (inner(grad(p), w) * dX) -
                             (inner(v, grad(q)) * dX);

         //Arif Masud and Thomas JR Hughes. A stabilized mixed finite element method for darcy flow
        auto stab = 0.5 * (
            inner(grad(p), grad(q)) * dX - inner(v, w)       * dX + 
            inner(v, grad(q))       * dX - inner(grad(p), w) * dX
        );

        USparseMatrix A;
        UVector rhs, x;

        x = local_zeros(V.dof_map().n_local_dofs());
        forcing_function.eval(x, rhs);

        double norm_rhs = norm2(rhs);
        std::cout << norm_rhs << std::endl;

        assemble(bilinear_form + stab, A);
        apply_boundary_conditions(V.dof_map(), A, rhs);

        Factorization<USparseMatrix, UVector>().solve(A, rhs, x);

        write("rhs.e", V, rhs);
        write("sol.e", V, x);

        auto v_h = interpolate(x, v);
        UVector divergence_v;

        assemble(inner(div(v_h), q) * dX, divergence_v);
        double measure_div = sum(divergence_v);
        std::cout <<  "measure_div: " << measure_div << std::endl;

        write("div.e", V, divergence_v);
    }
}
