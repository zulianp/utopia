#include "utopia_Intrepid2Test.hpp"
#include <algorithm>

#ifdef WITH_INTREPID2
#include <utopia_Intrepid2_Assembler.hpp>
#include "utopia_libmesh.hpp"

namespace utopia {

	void run_intrepid2_test(libMesh::LibMeshInit &init)
	{
		libMesh::DistributedMesh mesh(init.comm());

		int n = 10;
		libMesh::MeshTools::Generation::build_square(
			mesh,
			n,  n,
			0., 1.,
			0., 1.,
			libMesh::TRI3
		);

		LibMeshFunctionSpace V(mesh); //TODO

		auto u = trial(V);
		init_constraints( constraints(
			boundary_conditions(u == coeff(0), {1, 3}),
			boundary_conditions(u == coeff(0), {0, 2})
		));
		

		V.initialize();

		//////////////////////////////////////////
		//////////////////////////////////////////

		


		auto v = test(V);

		auto b_form = inner(grad(u), grad(v)) * dX;
		auto l_form = inner(coeff(1.), v) * dX;


		//USparseMatrix H;
		//UVector rhs;

		TSMatrixd H;
		TVectord rhs;

		Intrepid2Assembler assembler;
		assembler.assemble(l_form, rhs);
		assembler.assemble(b_form, H);

		disp(H);

		apply_boundary_conditions(V.dof_map(), H, rhs);

		TVectord x = local_zeros(local_size(rhs)); //UVector x
		solve(H, rhs, x);

		//write("intrepid2test.e", V, x); //TODO
		//////////////////////////////////////////
	}
}

#else
	
	namespace utopia {
		void run_intrepid2_test(libMesh::LibMeshInit &){}
	}

#endif //WITH_INTREPID2

