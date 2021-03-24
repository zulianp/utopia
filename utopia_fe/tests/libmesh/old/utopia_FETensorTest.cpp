#include "utopia_FETensorTest.hpp"

#include "utopia_SymbolicFunction.hpp"
#include "utopia_ui.hpp"

#include "utopia_Flow.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

#include "utopia_FEEval_Local.hpp"
#include "utopia_libmesh_AssembleLocal.hpp"

#include "libmesh/mesh_generation.h"
#include "utopia_libmesh.hpp"

#include "libmesh/boundary_mesh.h"
#include "libmesh/mesh_refinement.h"

namespace utopia {

    void FETensorTest::run(Input &in) {
        std::cout << "[FETensorTest]" << std::endl;
        using Space = utopia::LibMeshFunctionSpace;
        using FE = utopia::FiniteElement<Space>;

        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<Space> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        UIForcingFunction<Space, UVector> forcing_function(space.subspace(0));
        in.get("forcing-function", forcing_function);

        bool no_solve = false;
        in.get("no-solve", no_solve);

        auto &V = space.space().subspace(0);
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////

        // Global tensors
        USparseMatrix A;
        UVector rhs, x, forcing_term;

        // Local tensors
        FormVectord g;
        FormScalard f;
        FormVectord alpha_g;
        MultiScalard dx;

        std::cout << "ASSEMBLY begin" << std::endl;
        Chrono c;
        c.start();

        utopia::assemble(V,
                         // FIXME
                         inner(grad(trial(V)), grad(test(V))) * dX,
                         A,
                         [&](FE &element, USerialMatrix &mat) {
                             // Symbolic
                             auto u = trial(element);

                             dx = measure(element);
                             g = grad(u);
                             alpha_g = 0.5 * g;
                             mat = form2(alpha_g, g, dx);
                         });

        utopia::assemble(V,
                         // FIXME
                         inner(coeff(0.0), test(V)) * dX,
                         rhs,
                         [&](FE &element, USerialVector &vec) {
                             // Symbolic
                             auto u = trial(element);

                             // Symbolic to Numeric conversion
                             f = u;

                             f *= 0.0;
                             dx = measure(element);
                             vec = form1(f, dx);
                         });

        c.stop();

        std::cout << "ASSEMBLY end: " << c << std::endl;
        if (no_solve) return;

        x.zeros(row_layout(A));
        forcing_term.zeros(layout(x));
        forcing_function.eval(x, forcing_term);
        rhs += forcing_term;

        apply_boundary_conditions(V, A, rhs);

        Factorization<USparseMatrix, UVector> fact;
        fact.describe(std::cout);
        fact.solve(A, rhs, x);

        write("rhs.e", V, rhs);
        write("sol.e", V, x);
    }
}  // namespace utopia
