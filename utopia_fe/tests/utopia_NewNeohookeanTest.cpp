#include "utopia_NewNeohookeanTest.hpp"

#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_Flow.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_FEEval_MultiTensor.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"
#include "utopia_NewNeoHookean.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/mesh_generation.h"

#include "libmesh/mesh_refinement.h"
#include "libmesh/boundary_mesh.h"

namespace utopia {

    void NewNeohookeanTest::run(Input &in)
    {
        std::cout << "[NewNeohookeanTest]" << std::endl;
        using Space = utopia::LibMeshFunctionSpace;
        using FE = utopia::FiniteElement<Space>;
        using TFE = utopia::FiniteElement<ProductFunctionSpace<Space>>;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<Space> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<Space, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);

        auto &V = space.space();

        //FIXME
        double lambda = 1.0, mu = 1.0, rescaling = 1.0;

        //FIXME use ghost
        UVector x = local_zeros(V[0].dof_map().n_local_dofs());

        LameeParameters params;
        NewNeoHookean<ProductFunctionSpace<Space>, USparseMatrix, UVector> material(V, params);

        USparseMatrix H;
        UVector g;
        material.assemble_hessian_and_gradient(x, H, g);
    }

}