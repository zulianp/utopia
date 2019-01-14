#include "utopia_FractureFlowUtils.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_TransferAssembler.hpp"

#include "moonolith_communicator.hpp"

#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"

#include <vector>

namespace utopia {

	std::shared_ptr<SemiGeometricMultigrid> make_mg_solver(
		const LibMeshFunctionSpace &space, const int n_levels)
	{
		auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
		auto smoother      = std::make_shared<GaussSeidel<USparseMatrix, UVector>>();
		auto mg            = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);

		mg->algebraic().rtol(1e-9);
		mg->algebraic().atol(1e-14);
	    // mg->verbose(true);
		mg->init(space, n_levels);
		return mg;
	}

	void refine_around_fractures(
		const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
		const libMesh::Order &elem_order,
		const std::shared_ptr<libMesh::UnstructuredMesh> &mesh,
		const int refinement_loops,
		const bool use_interpolation
		)
	{

		libMesh::MeshRefinement mesh_refinement(*mesh);

		for(int i = 0; i < refinement_loops; ++i) {
	    //equations system
			auto vol_equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
			auto &vol_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("vol_sys");

			auto surf_equation_systems = std::make_shared<libMesh::EquationSystems>(*fracture_network);
			auto &surf_sys = surf_equation_systems->add_system<libMesh::LinearImplicitSystem>("surf_sys");

	    //scalar function space
			auto V_vol  = LibMeshFunctionSpace(vol_equation_systems, libMesh::LAGRANGE, elem_order,      "u_vol");
			auto V_surf = LibMeshFunctionSpace(surf_equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_surf");

			V_vol.initialize();
			V_surf.initialize();

			Chrono c;
			c.start();
			USparseMatrix B;
			moonolith::Communicator comm(mesh->comm().get());
			if(assemble_volume_transfer(
				comm,
				mesh,
				fracture_network,
				make_ref(V_vol.dof_map()),
				make_ref(V_surf.dof_map()),
				0,
				0,
				true,
				1,
				B,
				{},
				use_interpolation))
			{
				c.stop();
				std::cout << c << std::endl;

				Interpolator interp(make_ref(B));
				interp.normalize_rows();
				interp.describe(std::cout);

				USparseMatrix D_inv = diag(1./sum(B, 1));
				USparseMatrix T = D_inv * B;

				USparseMatrix T_t = transpose(T);
				UVector t_temp = sum(T_t, 1);
				UVector t = ghosted(local_size(t_temp).get(0), size(t_temp).get(0), V_vol.dof_map().get_send_list());
				t = t_temp;


				std::vector<libMesh::dof_id_type> indices;
				std::vector<double> values;

				mesh_refinement.clean_refinement_flags();

				Read<UVector> r_(t);
				for(auto e_it = elements_begin(*mesh); e_it != elements_end(*mesh); ++e_it) {
					V_vol.dof_map().dof_indices(*e_it, indices);
					t.get(indices, values);

					double val = std::accumulate(values.begin(), values.end(), 0.,  std::plus<double>());
					if(val > 0) {
						(*e_it)->set_refinement_flag(libMesh::Elem::REFINE);
					}

				}

				mesh_refinement.make_flags_parallel_consistent();
				mesh_refinement.refine_elements();
				mesh_refinement.test_level_one(true);

			} else {
				assert(false);
			}
		}

	    // mesh_refinement.clean_refinement_flags();
		mesh->prepare_for_use();
	}

	void write_solution(
		const std::string &name,
		UVector &sol,
		const int time_step,
		const double t,
		LibMeshFunctionSpace &space,
		libMesh::Nemesis_IO &io)
	{
		utopia::convert(sol, *space.equation_system().solution);
		space.equation_system().solution->close();
		io.write_timestep(name, space.equation_systems(), time_step, t);
	}

	void transform_values(
		const LibMeshFunctionSpace &from,
		const UVector &from_vector,
		const LibMeshFunctionSpace &to, 
		UVector &to_vector,
		std::function<double(const double &)> fun)
	{
		using libMesh::dof_id_type;


		if(empty(to_vector)) {
			to_vector = local_zeros(to.dof_map().n_local_dofs());
		}

		const auto &m = from.mesh();
		assert(&m == &to.mesh());

		std::vector<dof_id_type> from_dofs, to_dofs;
		
		Read<UVector> r_(from_vector);
		Write<UVector> w_(to_vector);

		for(auto e_it = elements_begin(m); e_it != elements_end(m); ++e_it) {
			const auto &e = *e_it;

			from.dof_map().dof_indices(e, from_dofs, from.subspace_id());
			to.dof_map().dof_indices(e, to_dofs, to.subspace_id());

			assert(from_dofs.size() == to_dofs.size());

			auto n = from_dofs.size();

			for(std::size_t i = 0; i < n; ++i) {
				to_vector.set(to_dofs[i], fun(from_vector.get(from_dofs[i])));
			}
		}
	}


	void copy_values(
		const LibMeshFunctionSpace &from,
		const UVector &from_vector,
		const LibMeshFunctionSpace &to, 
		UVector &to_vector)
	{

		transform_values(from, from_vector, to, to_vector, [](const double &value) -> double { return value; });
		// using libMesh::dof_id_type;


		// if(empty(to_vector)) {
		// 	to_vector = local_zeros(to.dof_map().n_local_dofs());
		// }

		// const auto &m = from.mesh();
		// assert(&m == &to.mesh());

		// std::vector<dof_id_type> from_dofs, to_dofs;
		
		// Read<UVector> r_(from_vector);
		// Write<UVector> w_(to_vector);

		// for(auto e_it = elements_begin(m); e_it != elements_end(m); ++e_it) {
		// 	const auto &e = *e_it;

		// 	from.dof_map().dof_indices(e, from_dofs, from.subspace_id());
		// 	to.dof_map().dof_indices(e, to_dofs, to.subspace_id());

		// 	assert(from_dofs.size() == to_dofs.size());

		// 	auto n = from_dofs.size();

		// 	for(std::size_t i = 0; i < n; ++i) {
		// 		to_vector.set(to_dofs[i], from_vector.get(from_dofs[i]));
		// 	}
		// }
	}

}
