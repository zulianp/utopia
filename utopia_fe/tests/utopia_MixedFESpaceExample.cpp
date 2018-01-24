#include "utopia_MixedFESpaceExample.hpp"

#include <iostream>
#include "utopia.hpp"

//fe extension
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>

#include "utopia_ContactSimParams.hpp"
#include "utopia_MixedFESpace.hpp"

using namespace utopia;
using namespace std;
using namespace libMesh;


namespace utopia {
	void run_mixed_fe_space_example(libMesh::LibMeshInit &init)
	{
	// 	auto mesh = make_shared<Mesh>(init.comm());		
	// 	MeshTools::Generation::build_cube (*mesh,
	// 		10, 10, 10,
	// 		-1., 1.,
	// 		-1., 1.,
	// 		-1., 1.,
	// 		HEX27);

	// 	const int dim = mesh->mesh_dimension();

	// 	LibMeshFEContext<LinearImplicitSystem> context(mesh);
	// 	auto Uh1 = fe_space(LAGRANGE, FIRST, context);
	// 	auto Uh2 = fe_space(LAGRANGE, SECOND, context);
	// 	auto Uh3 = fe_space(LAGRANGE, THIRD, context);

	// 	auto mixed_space = Uh1 * Uh2 * Uh3;
	// 	mixed_space.each(PrintIndex());

	// 	auto U1 = mixed_space.sub<0>();
	// 	auto U2 = mixed_space.sub<1>();
	// 	auto U3 = mixed_space.sub<2>();

	// 	auto u1    = fe_function(U1);
	// 	auto m_fun = fe_function(mixed_space);

	// 	auto bf = integral(dot(m_fun, m_fun));

	// 	std::cout << tree_format(bf.getClass()) << std::endl;
	}
}