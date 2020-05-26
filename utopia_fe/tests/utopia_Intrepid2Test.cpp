#include "utopia_Intrepid2Test.hpp"
#include <algorithm>

#include "utopia_fe_kokkos_fix.hpp"

#include "utopia_Base.hpp"
#include "utopia_fe_base.hpp"

// TODO replace them
//#include "utopia_libmesh_FunctionSpace.hpp"
//#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEFunction.hpp"

#ifdef WITH_INTREPID2
#include <utopia_Intrepid2_Assembler.hpp>
#include "utopia_Intrepid2Backend.hpp"

#include <libmesh/mesh_generation.h>
#include "libmesh/parallel_mesh.h"

namespace utopia {
    void Intrepid2Test::run(Input &in) {
        // FIXME add a real communicator
        int argc;
        const char *const *argv;
        libMesh::LibMeshInit init(argc, argv);
        libMesh::DistributedMesh mesh(init.comm());
        //
        int n = 10;
        libMesh::MeshTools::Generation::build_square(mesh, n, n, 0., 1., 0., 1., libMesh::TRI3);
        std::string var_name = "pota";
        //		Intrepid2FunctionSpace V(var_name);//mesh); //TODO

        //		auto u = trial(V);
        /*init_constraints( constraints(
                boundary_conditions(u == coeff(0), {1, 3}),
                boundary_conditions(u == coeff(0), {0, 2})
        ));*/

        //		V.initialize();

        //////////////////////////////////////////
        //////////////////////////////////////////
        //		auto v = test(V);

        //		auto b_form = inner(grad(u), grad(v)) * dX;
        //		auto l_form = inner(coeff(1.), v) * dX;

        TSUSerialMatrix H;
        TUSerialVector rhs;

        Intrepid2Assembler assembler;
        //		assembler.assemble(l_form, rhs);
        //		assembler.assemble(b_form, H);

        disp(H);
        // TODO FIXME
        //		apply_boundary_conditions(V.dof_map(), H, rhs);

        TUSerialVector x = local_zeros(local_size(rhs));  // UVector x
        solve(H, rhs, x);

        //		write("intrepid2test.e", V, x); //TODO implement write function
        //////////////////////////////////////////
    }
}  // namespace utopia

#else

namespace utopia {
    void Intrepid2Test::run(Input &in) {}
}  // namespace utopia

#endif  // WITH_INTREPID2
