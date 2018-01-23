#include "utopia_SemiGeometricMultigrid.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_Socket.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "moonolith_communicator.hpp"


#include "libmesh/mesh_tools.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/libmesh_version.h"

#include <cmath>

namespace utopia {
	SemiGeometricMultigrid::SemiGeometricMultigrid()
	: mg(
		std::make_shared<GaussSeidel<DSMatrixd, DVectord>>(),
		std::make_shared<Factorization<DSMatrixd, DVectord>>()
		)
	{ }

	void SemiGeometricMultigrid::init(const libMesh::EquationSystems &es, const std::size_t n_levels)
	{
		const int system_number = 0;
		const auto &mesh = es.get_mesh();
		const auto &dof_map = es.get_system(system_number).get_dof_map();
		const auto dim = mesh.mesh_dimension();

		mg.set_fix_semidefinite_operators(true);

#if LIBMESH_VERSION_LESS_THAN(1, 0, 3)
		libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::bounding_box(mesh);
#else
		libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::create_bounding_box(mesh);
#endif

		auto r = bb.max() - bb.min();

		const std::size_t n_coarse_spaces = n_levels - 1;
		meshes.resize(n_coarse_spaces);
		equation_systems.resize(n_coarse_spaces);

		switch(dim) {
			case 2: 
			{
				const int n_segments = std::max(2, int( ceil( std::sqrt(mesh.n_active_elem() / std::pow(4, n_coarse_spaces)) ) ) );
				const double aspect_ratio = r(0)/r(1);
				const int nx = std::max(int(ceil(n_segments * aspect_ratio)), 2);
				const int ny = std::max(n_segments, 2);

				auto m = std::make_shared<libMesh::DistributedMesh>(mesh.comm());
				libMesh::MeshTools::Generation::build_square (
					*m,
					nx, ny,
					bb.min()(0), bb.max()(0),
					bb.min()(1), bb.max()(1),
					libMesh::QUAD4);

				meshes[0] = m;
				equation_systems[0] = std::make_shared<libMesh::EquationSystems>(*m);

				// plot_mesh(*m, "mg/l_" + std::to_string(0));

				//use refinement instead
				for(std::size_t i = 1; i < n_levels-1; ++i) {
					auto m_i = std::make_shared<libMesh::DistributedMesh>(*meshes[i-1]);
					
					{
						libMesh::MeshRefinement mesh_refinement(*m_i);
						mesh_refinement.make_flags_parallel_consistent();
						mesh_refinement.uniformly_refine(1);
					}

					equation_systems[i] = std::make_shared<libMesh::EquationSystems>(*m_i);
					meshes[i] = m_i;

					// plot_mesh(*m_i, "mg/l_" + std::to_string(i));
				}

				break;
			}

			default:
			{	
				assert(false && "implement me");
				break;
			}
		}

		// plot_mesh(mesh, "mg/L");

		for(std::size_t i = 0; i < n_levels-1; ++i) {
			auto &sys = equation_systems[i]->add_system<libMesh::LinearImplicitSystem>(es.get_system(system_number).name());

			for(unsigned int v = 0; v < dof_map.n_variables(); ++v) {
				const auto &var = dof_map.variable(v);
				//FIXME
				sys.add_variable(var.name(), dof_map.variable_order(v), libMesh::LAGRANGE);
			}

			sys.init();
		}

		std::vector<std::shared_ptr<DSMatrixd>> interpolators(n_coarse_spaces);
		moonolith::Communicator comm(mesh.comm().get());
		
		for(std::size_t i = 1; i < n_coarse_spaces; ++i) {
			interpolators[i-1] = std::make_shared<DSMatrixd>();

			bool success = assemble_volume_transfer(
				comm,
				meshes[i-1],
				meshes[i],
				make_ref(equation_systems[i-1]->get_system(0).get_dof_map()), 
				make_ref(equation_systems[i]->get_system(0).get_dof_map()),
				0,
				0,
				true,
				dof_map.n_variables(),
				*interpolators[i-1]
				); assert(success);
		}

		interpolators[n_coarse_spaces-1] = std::make_shared<DSMatrixd>();
		bool success = assemble_volume_transfer(
			comm,
			meshes[n_coarse_spaces-1],
			make_ref(const_cast<libMesh::MeshBase &>(mesh)),
			make_ref(equation_systems[n_coarse_spaces-1]->get_system(0).get_dof_map()), 
			make_ref(const_cast<libMesh::DofMap &>(dof_map)),
			0,
			0,
			true,
			dof_map.n_variables(),
			*interpolators[n_coarse_spaces-1]
			); assert(success);

		if(mg.verbose()) {
			for(const auto &e : equation_systems) {
				std::cout << "dofs: " << e->get_system(0).get_dof_map().n_dofs() << std::endl;
			}

			std::cout << "dofs: " << es.get_system(0).get_dof_map().n_dofs() << std::endl;
		}

		mg.init_transfer_from_fine_to_coarse(std::move(interpolators));
		//FIXME naming is wrong
		// mg.init_transfer_from_coarse_to_fine(std::move(interpolators));
	}	

	void SemiGeometricMultigrid::update(const std::shared_ptr<const DSMatrixd> &op)
	{
		mg.galerkin_assembly(op);
	}

	bool SemiGeometricMultigrid::apply(const DVectord &rhs, DVectord &sol)
	{
		return mg.solve(rhs, sol);
	}

}