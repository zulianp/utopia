#include "utopia_SemiGeometricMultigrid.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_Socket.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_Contact.hpp"

#include "moonolith_communicator.hpp"


#include "libmesh/mesh_tools.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/libmesh_version.h"

#include <cmath>

namespace utopia {
	static void make_d(const DSMatrixd &mat, DVectord &res)
	{
		res = sum(mat, 1);

		ReadAndWrite<DVectord> rw_(res);
		auto r = range(res);
		for(auto k = r.begin(); k != r.end(); ++k) {
			if(approxeq(res.get(k), 0.0, 1e-14)) {
				res.set(k, 1.);
			}
		}
	}

	SemiGeometricMultigrid::SemiGeometricMultigrid(
		const std::shared_ptr<Smoother<DSMatrixd, DVectord> > &smoother,
		const std::shared_ptr<LinearSolver<DSMatrixd, DVectord> > &linear_solver)
	: mg(smoother, linear_solver), is_block_solver_(false)
	{ }

	void SemiGeometricMultigrid::init(const libMesh::EquationSystems &es, const std::size_t n_levels)
	{
		const int system_number = 0;
		const auto &mesh = es.get_mesh();
		const auto &dof_map = es.get_system(system_number).get_dof_map();
		const auto dim = mesh.mesh_dimension();

		mg.set_fix_semidefinite_operators(true);

#if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
		libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::bounding_box(mesh);
#else
		libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::create_bounding_box(mesh);
#endif

		auto r = bb.max() - bb.min();

		const std::size_t n_coarse_spaces = n_levels - 1;
		meshes.resize(n_coarse_spaces);
		equation_systems.resize(n_coarse_spaces);

		auto m = std::make_shared<libMesh::DistributedMesh>(mesh.comm());

		switch(dim) {
			case 2: 
			{
				const int n_segments = std::max(2, int( std::round( std::sqrt(mesh.n_nodes() / std::pow(4, n_coarse_spaces)) ) ) );
				const double max_r = std::max(r(0), r(1));

				const double aspect_ratio_x = r(0)/max_r;
				const double aspect_ratio_y = r(1)/max_r;

				const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
				const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);

				libMesh::MeshTools::Generation::build_square (
					*m,
					nx, ny,
					bb.min()(0), bb.max()(0),
					bb.min()(1), bb.max()(1),
					libMesh::QUAD4);
				
				break;
			}

			case 3:
			{
				const int n_segments = std::max(2, int( std::round( std::cbrt(mesh.n_nodes() / std::pow(8, n_coarse_spaces)) ) ) );
				
				const double max_r = std::max(r(0), std::max(r(1), r(2)));
				
				const double aspect_ratio_x = r(0)/max_r;
				const double aspect_ratio_y = r(1)/max_r;
				const double aspect_ratio_z = r(2)/max_r;

				const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
				const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);
				const int nz = std::max(int(std::round(n_segments * aspect_ratio_z)), 2);

				libMesh::MeshTools::Generation::build_cube (
					*m,
					nx, ny, nz,
					bb.min()(0), bb.max()(0),
					bb.min()(1), bb.max()(1),
					bb.min()(2), bb.max()(2),
					libMesh::HEX8);

				break;
			} 

			default:
			{	
				assert(false && "implement me");
				break;
			}
		}

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

		// plot_mesh(mesh, "mg/L");

		for(std::size_t i = 0; i < n_levels-1; ++i) {
			auto &sys = equation_systems[i]->add_system<libMesh::LinearImplicitSystem>(es.get_system(system_number).name());

			for(unsigned int v = 0; v < dof_map.n_variables(); ++v) {
				const auto &var = dof_map.variable(v);
				//FIXME
				sys.add_variable(var.name(), libMesh::FIRST, libMesh::LAGRANGE);
			}

			sys.init();
		}

		interpolators_.resize(n_coarse_spaces);
		moonolith::Communicator comm(mesh.comm().get());

		DVectord d_diag;
		
		for(std::size_t i = 1; i < n_coarse_spaces; ++i) {
			interpolators_[i-1] = std::make_shared<DSMatrixd>();

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
				*interpolators_[i-1]
				); assert(success);

			DVectord d_diag;

			make_d(*interpolators_[i-1], d_diag);
			*interpolators_[i-1] = diag(1./d_diag) * *interpolators_[i-1];

		}

		interpolators_[n_coarse_spaces-1] = std::make_shared<DSMatrixd>();
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
			*interpolators_[n_coarse_spaces-1]
			); assert(success);


		make_d(*interpolators_[n_coarse_spaces-1], d_diag);
		*interpolators_[n_coarse_spaces-1] = diag(1./d_diag) * *interpolators_[n_coarse_spaces-1];

		// write("T.m", *interpolators[n_coarse_spaces-1]);

		if(mg.verbose()) {
			for(const auto &e : equation_systems) {
				std::cout << "dofs: " << e->get_system(0).get_dof_map().n_dofs() << std::endl;
			}

			std::cout << "dofs: " << es.get_system(0).get_dof_map().n_dofs() << std::endl;
		}

		mg.init_transfer_from_fine_to_coarse(interpolators_);
		//FIXME naming is wrong
		// mg.init_transfer_from_coarse_to_fine(std::move(interpolators));
	}	

	void SemiGeometricMultigrid::update(const std::shared_ptr<const DSMatrixd> &op)
	{
		mg.update(op);

		//hacky
		if(is_block_solver_) {
			for(SizeType i = 0; i < mg.num_levels(); ++i) {
				const_cast<DSMatrixd &>(mg.level(i).A()).implementation().convert_to_mat_baij(meshes[0]->mesh_dimension());
			}
		}

	}

	bool SemiGeometricMultigrid::apply(const DVectord &rhs, DVectord &sol)
	{
		return mg.apply(rhs, sol);
	}

	void SemiGeometricMultigrid::update_contact(Contact &contact)
	{
		const auto last_interp = mg.num_levels() - 2;
		auto c_I = std::make_shared<DSMatrixd>();
		*c_I = transpose(contact.complete_transformation) * *interpolators_[last_interp];
		mg.update_transfer(last_interp, Transfer<DSMatrixd, DVectord>(c_I));
	}

}
