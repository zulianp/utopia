#include "utopia_SDCTest.hpp"

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

#include "libmesh/exodusII_io.h"
#include <algorithm>


#include <memory>


namespace utopia {

	template<class Matrix, class Vector>
	class ExplicitSDC {
	public:
		ExplicitSDC()
		: max_iter_(10)
		{}

		template<class Fun>
		void integrate(const std::vector<double> &quad_points,
					   const std::vector<double> &quad_weights,
					   const double tbegin,
					   const double tend,
					   const int n_time_steps,
					   const Fun &fun,
					   const Matrix &mass_matrix,
					   const Vector &u0,
					   std::vector<Vector> &solutions
					   )
		{
			auto s = local_size(u0);

			const std::size_t N = quad_points.size();
			const double dt = (tend - tbegin)/n_time_steps;

			solutions.resize(N);
			corrections_.resize(N);


			for(auto &s : solutions) {
				s = local_zeros(s);
			}

			for(auto &c : corrections_) {
				c = local_zeros(s);
			}


			for(std::size_t i = 0; i < max_iter_; ++i) {

				for(std::size_t k = 1; k < N; ++k) {
					const double dtk  = quad_points[k] - quad_points[k-1];
					const double wk   = quad_weights[k];
					const double wkm1 = quad_weights[k-1];
				}
			}
		}

		std::size_t max_iter_;
		std::vector<Vector> corrections_;

	};

	void run_sdc_test(libMesh::LibMeshInit &init)
	{
		auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());

		const unsigned int n = 10;
		libMesh::MeshTools::Generation::build_square(*lm_mesh,
			n, n,
			0, 1,
			0, 1.,
			libMesh::QUAD8);

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*lm_mesh);
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("sdc");


		auto V = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u");
		auto u = trial(V);
		auto v = test(V);


		const double dt = 0.01;
		const std::size_t n_ts = 40;

		UVector storage_u;
		auto R = coeff(1.);
		auto var_u = interpolate(storage_u, u);
		auto form =
				( inner(R, v) - inner( grad(var_u), grad(v)) ) * dX;

		auto constr = constraints(
						boundary_conditions(u == coeff(0.),  {1, 3}),
						boundary_conditions(u == coeff(0.),  {0}),
						boundary_conditions(u == coeff(0.0), {2})
					);

		//FIXME
		FEBackend<LIBMESH_TAG>::init_constraints(constr);
		V.initialize();

		UVector u0 = local_zeros(V.dof_map().n_local_dofs());
		//apply boundary conditions

		auto f = [&form, &storage_u](const UVector &u, UVector &ret) {
			storage_u = u;
			// assemble(form, ret);
		};

		std::vector<double> quad_points{0, 0.5, 1.};
		ExplicitSDC<USMatrix, UVector> sdc;
	}
}
