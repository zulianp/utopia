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
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));

        in.get("mesh", mesh);
        in.get("space", space);

        auto dim = mesh.mesh().mesh_dimension();
        auto &M  = space.space();
        
        if(dim + 1 != M.n_subspaces()) {
            std::cerr << "[Error] wrong number of subspaces, expected " <<   (dim + 1) << std::endl;
            std::cerr << "V[0] for scalar, V[1, ..., dim + 1] for gradient" << std::endl;
            return;
        }

        auto &V = M.subspace(0);
        auto W  = M.subspace(1, dim + 1);

        auto u = trial(V);
        auto v = test(V);

        auto sigma = trial(W);
        auto tau   = test(W);
        
        auto tau_x = test(W[0]);
        auto tau_y = test(W[1]);

        //
        auto a   = -(inner(sigma, tau) * dX);
        auto b_t = inner(u, div(tau)) * dX;

        auto b   = inner(div(sigma), v) * dX;


        auto bilinear_form = a + b + b_t;

        auto fv   = inner(coeff(0.0), v) * dX;

        auto ftau = inner(coeff(0.0), tau_x) * dX + inner(coeff(0.0), tau_y) * dX;

        auto linear_form = fv + ftau;

        USparseMatrix A;
        UVector rhs, x;

        assemble(bilinear_form == linear_form, A, rhs);

        write("A_neu.m", A);

        apply_boundary_conditions(V.dof_map(), A, rhs);

        write("A.m", A);

        x = local_zeros(local_size(rhs));

        Factorization<USparseMatrix, UVector>().solve(A, rhs, x);

        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}
