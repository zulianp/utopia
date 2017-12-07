#include "utopia_MechTest.hpp"

#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_Mechanics.hpp"

#include "libmesh/exodusII_io.h"
#include <algorithm>


#include <memory>

namespace utopia {

	void run_mech_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		
		const unsigned int n = 8;
		const unsigned int dim = 2;
		const double dt = 1;

		libMesh::MeshTools::Generation::build_square(
			*mesh,
			n, n,
			0, 1,
			0, 1.,
			libMesh::QUAD8);

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("mech_test");

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::SECOND, "disp_y");
		auto V = Vx * Vy;

		auto ux = trial(Vx);
		auto uy = trial(Vy);

		auto constr = constraints(
			boundary_conditions(ux == coeff(0.),  {0, 2}),
			boundary_conditions(uy == coeff(.1), {0}),
			boundary_conditions(uy == coeff(-.1),  {2})
		);

		FEBackend<LIBMESH_TAG>::init_constraints(constr);

		Vx.initialize();

		DVectord u = local_zeros(Vx.dof_map().n_local_dofs());
		DVectord internal_force;
		DSMatrixd stiffness_matrix;

		MechanicsContext mech_ctx;
		MechanicsState old;
		MechanicsState current;

		mech_ctx.init_mass_matrix(V);

		old.init(local_size(u), size(u));
		current.init(local_size(u), size(u));

		auto elast = std::make_shared<LinearElasticity>();

		LameeParameters params;
		elast->init(V, params, u, mech_ctx.stiffness_matrix, old.internal_force);

		auto v = test(V);
		auto vx = test(Vx);
		auto vy = test(Vy);
		
		auto ef = std::make_shared<ConstantExternalForce>();

		LMDenseVector f_value = values(mesh->mesh_dimension(), 0.);

		// ef->init(inner(coeff(f_value), v) * dX);
		ef->init((inner(coeff(0.5), vx) + inner(coeff(-0.1), vy)) * dX);

		apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, ef->value);
		ef->update(old.t, old.external_force);
		ef->update(old.t + dt, current.external_force);

		disp(mech_ctx.stiffness_matrix);
		disp(old.external_force);

		auto integrator = std::make_shared<ImplicitEuler>(dim, Vx.dof_map());
		integrator->apply(dt, mech_ctx, old, current);

		disp(current.displacement);

		convert(current.displacement, *sys.solution);
		sys.solution->close();

		libMesh::ExodusII_IO(*mesh).write_equation_systems ("mech_test.e", *equation_systems);
	}
}
