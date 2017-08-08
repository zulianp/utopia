#include "utopia_GeometryTest.hpp"

#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include "utopia_fe.hpp"
#include "utopia_LibMeshBackend.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>

#include <memory>

namespace utopia {
	void run_geometry_test(libMesh::LibMeshInit &init)
	{
		using namespace libMesh;
		using namespace std;

		const auto order_elem = FIRST;

		auto mesh = make_shared<Mesh>(init.comm());
		mesh->read("/Users/patrick/Desktop/PostDOC/sccer_turbines/turbine.e");

		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		auto space_1 = fe_space(LAGRANGE, order_elem, context);
		auto space_2 = fe_space(LAGRANGE, order_elem, context);
		auto space_3 = fe_space(LAGRANGE, order_elem, context);
		
		context.equation_systems.init();

		DVectord is_normal_component;
		DVectord normals;
		DSMatrixd mat;
		assemble_normal_tangential_transformation(*mesh, space_1.dof_map(), {5}, is_normal_component, normals, mat);
	}
}