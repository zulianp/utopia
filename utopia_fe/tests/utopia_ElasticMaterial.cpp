#include "utopia_ElasticMaterial.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"

#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

// #include "utopia_Equations.hpp"
// #include "utopia_FEConstraints.hpp"
// #include "utopia_FindSpace.hpp"
// #include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
// #include "utopia_FEKernel.hpp"

// #include "libmesh/exodusII_io.h"
#include <algorithm>


namespace utopia {
	typedef utopia::ProductFunctionSpace<LibMeshFunctionSpace> FunctionSpaceT;

	bool LinearElasticity::init(
		const FunctionSpaceT &V,
		const LameeParameters &params,
		const DVectord &displacement0,
		DSMatrixd &stiffness_matrix,
		DVectord  &internal_stress)
	{
		auto u = trial(V);
		auto v = test(V);

		auto mu     = params.var_mu();
		auto lambda = params.var_lambda();

		auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
		auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

		auto b_form = integral((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v)));

		if(!assemble(b_form, stiffness_matrix)) return false;
		return update(displacement0, stiffness_matrix, internal_stress);
	}

	bool LinearElasticity::update(
			const DVectord &displacement,
			DSMatrixd &stiffness_matrix,
			DVectord  &internal_stress)
	{
		internal_stress = stiffness_matrix * displacement;
		return true;
	}
}
