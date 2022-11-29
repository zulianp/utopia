#include "utopia_MassApp.hpp"

#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_ui.hpp"

#include "utopia_MeshTransferOperator.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh_refinement.h"

namespace utopia {

    void MassApp::run(Input &in) {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));

        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<LibMeshFunctionSpace, UVector> forcing_function(space.subspace(0));

        in.get("forcing-function", forcing_function);

        auto &V = space.space().subspace(0);

        auto u = trial(V);
        auto v = test(V);

        Adaptivity a;

        USparseMatrix pre_constraint;

        USparseMatrix post_constraint;

        a.constraint_matrix(V.mesh(), V.dof_map(), pre_constraint, post_constraint);

        USparseMatrix Mass;
        UVector rhs, x, forcing_term;

        auto linear_form = inner(coeff(0.0), v) * dX;

        auto bilinear_form = inner(u, v) * dX;

        rhs = local_zeros(V.dof_map().n_local_dofs());

        x.zeros(layout(rhs));

        utopia::assemble(inner(u, v) * dX, Mass);

        forcing_term.zeros(layout(rhs));

        forcing_function.eval(x, forcing_term);

        // utopia::disp(forcing_term);

        rhs += forcing_term;

        rhs += post_constraint * rhs;

        Factorization<USparseMatrix, UVector> fact(MATSOLVERMUMPS, PCLU);
        fact.describe(std::cout);
        fact.solve(Mass, rhs, x);

        UVector sum_row = sum(Mass, 1);

        UVector m_inv = 1. / sum_row;

        UVector test = e_mul(m_inv, forcing_term);

        utopia::write("post.m", post_constraint);

        utopia::write("pre.m", pre_constraint);

        test += post_constraint * test;

        write("rhs.e", V, rhs);

        write("test.e", V, test);

        write("rhs.e", V, rhs);

        write("sol.e", V, x);
    }
}  // namespace utopia
