#include "utopia_TestVolume2SurfaceTransfer.hpp"

#include "utopia_libmesh.hpp"
#include "moonolith_communicator.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_TransferAssembler.hpp"

#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/mesh_refinement.h"
#include <algorithm>


typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

namespace utopia {
	void refine_around_fractures(
		const std::shared_ptr<libMesh::UnstructuredMesh> &fracture_network,
		const libMesh::Order &elem_order,
		const std::shared_ptr<libMesh::UnstructuredMesh> &mesh,
		const int refinement_loops = 1
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
			auto V_vol  = FunctionSpaceT(vol_equation_systems, libMesh::LAGRANGE, elem_order,      "u_vol");
			auto V_surf = FunctionSpaceT(surf_equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_surf");

			V_vol.initialize();
			V_surf.initialize();

			Chrono c;
			c.start();
			DSMatrixd B;
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
				true))
			{
				c.stop();
				std::cout << c << std::endl;

				Interpolator interp(make_ref(B));
				interp.normalize_rows();
				interp.describe(std::cout);

				DSMatrixd D_inv = diag(1./sum(B, 1));
				DSMatrixd T = D_inv * B;

				DSMatrixd T_t = transpose(T);
				DVectord t_temp = sum(T_t, 1);
				DVectord t = ghosted(local_size(t_temp).get(0), size(t_temp).get(0), V_vol.dof_map().get_send_list());
				t = t_temp;


				std::vector<libMesh::dof_id_type> indices;
				std::vector<double> values;

				mesh_refinement.clean_refinement_flags();

				Read<DVectord> r_(t);
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


	void run_volume_to_surface_transfer_test(libMesh::LibMeshInit &init)
	{
		auto n = 10;
		// auto elem_type  = libMesh::TET10;
		auto elem_type  = libMesh::TET4;
		// auto elem_type  = libMesh::HEX8;
		
		auto elem_order = libMesh::FIRST;
		// auto elem_order = libMesh::SECOND;

		// bool is_test_case = true;
		bool is_test_case = false;

		auto vol_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	
		auto surf_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());	

		// auto vol_mesh = std::make_shared<libMesh::ReplicatedMesh>(init.comm());	
		// auto surf_mesh = std::make_shared<libMesh::ReplicatedMesh>(init.comm());	


		if(is_test_case) {
			libMesh::MeshTools::Generation::build_cube(
				*vol_mesh,
				n, n, n,
				// 0., 1.,
				// 0., 1.,
				// -0.5, 0.5,
				-0.1, 1.1,
				-0.1, 1.1,
				-0.6, 0.6,
				elem_type
				);
			// surf_mesh->read("../data/test/square_with_2_tri.e");
			surf_mesh->read("../data/test/simple_network.e");
			// surf_mesh->set_mesh_dimension(3);


			libMesh::MeshRefinement mesh_refinement(*surf_mesh);
			mesh_refinement.make_flags_parallel_consistent();
			mesh_refinement.uniformly_refine(3);

			// refine_around_fractures(surf_mesh, elem_order, vol_mesh, 3);
		} else {
			
			// libMesh::MeshTools::Generation::build_cube(
			// 	*vol_mesh,
			// 	n, n, n,
			// 	-7.6, 7.6,
			// 	-7.6, 7.6,
			// 	-7.6, 7.6,
			// 	elem_type
			// 	);

			// surf_mesh->read("../data/test/fractures.e");


			libMesh::MeshTools::Generation::build_square(
				*vol_mesh,
				n, n,
				-0.5, 0.5,
				-0.5, 0.5,
				libMesh::TRI3
				);

			// vol_mesh->read("../data/frac/frac1d_background.e");
			surf_mesh->read("../data/frac/frac1d_network.e");
			
			{
				libMesh::MeshRefinement mesh_refinement(*surf_mesh);
				mesh_refinement.make_flags_parallel_consistent();
				mesh_refinement.uniformly_refine(4);
			}

			{
				refine_around_fractures(surf_mesh, elem_order, vol_mesh, 5);

				// libMesh::MeshRefinement mesh_refinement(*vol_mesh);
				// mesh_refinement.make_flags_parallel_consistent();
				// mesh_refinement.uniformly_refine(6);
			}
		}

		//equations system
		auto vol_equation_systems = std::make_shared<libMesh::EquationSystems>(*vol_mesh);
		auto &vol_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("vol_sys");
		auto &aux_sys = vol_equation_systems->add_system<libMesh::LinearImplicitSystem>("aux_sys");

		auto surf_equation_systems = std::make_shared<libMesh::EquationSystems>(*surf_mesh);
		auto &surf_sys = surf_equation_systems->add_system<libMesh::LinearImplicitSystem>("surf_sys");

		//scalar function space
		auto V_vol  = FunctionSpaceT(vol_equation_systems, libMesh::LAGRANGE, elem_order, "u_vol");
		auto V_aux  = FunctionSpaceT(vol_equation_systems, libMesh::LAGRANGE, elem_order, "indicator", aux_sys.number());

		auto V_surf = FunctionSpaceT(surf_equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_surf");

		V_vol.initialize();
		V_surf.initialize();
		V_aux.initialize();


		Chrono c;
		c.start();

		auto B = std::make_shared<DSMatrixd>();
		moonolith::Communicator comm(init.comm().get());

		const bool use_interpolation = true;
		if(assemble_volume_transfer(
			comm,
			vol_mesh,
			surf_mesh,
			make_ref(V_vol.dof_map()),
			make_ref(V_surf.dof_map()),
			0,
			0,
			true, 
			1,
			*B,
			{},
			use_interpolation))
		{

			c.stop();
			std::cout << c << std::endl;

			std::shared_ptr<TransferOperator> T;

			if(use_interpolation) {
				auto p_T = std::make_shared<Interpolator>(B);
				//only necessary for non-conforming mesh in parallel
				p_T->normalize_rows();
				T = p_T;
			} else {

				auto p_T = std::make_shared<PseudoL2TransferOperator>();
				p_T->init_from_coupling_operator(*B);
				T = p_T;

				DSMatrixd surf_mass_mat;
				assemble(inner(trial(V_surf), test(V_surf)) * dX, surf_mass_mat);
				double expected_mass = sum(surf_mass_mat);
				double actual_mass = sum(*B);

				std::cout << "operator volume : " << expected_mass << " == " <<  actual_mass << std::endl; 
			}

			T->describe(std::cout);

			DVectord v_vol = local_values(V_vol.dof_map().n_local_dofs(), 1.);
			// {
			// 	Write<DVectord> w_(v_vol);
			// 	v_vol.set(0, 0.);
			// }

			auto f_rhs = ctx_fun< std::vector<double> >([](const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> { 
				const auto &pts = ctx.fe()[0]->get_xyz();

				const auto n = pts.size();
				std::vector<double> ret(n);

				for(std::size_t i = 0; i != n; ++i) {
					double x = pts[i](0) - 0.5;
					double y = pts[i](1) - 0.5;
					double z = pts[i](2);
					ret[i] = std::abs(sin(x))*(x*x + y*y + z*z);
				}

				return ret;
			});

			auto u = trial(V_vol);
			auto v = test(V_vol);
			auto p_form = inner(f_rhs, v) * dX;
			auto m_form = inner(u, v) * dX;

			DVectord scaled_sol;
			DSMatrixd mass_mat;

			utopia::assemble(p_form, scaled_sol);
			utopia::assemble(m_form, mass_mat);

			if(elem_order == libMesh::FIRST) {
				v_vol = e_mul(1./sum(mass_mat, 1), scaled_sol);
			} else {
				Factorization<DSMatrixd, DVectord>().solve(mass_mat, scaled_sol, v_vol);
			}

			// v_vol.set(1.);

			double max_master = max(v_vol);

			v_vol *= 1./max_master;
			max_master = 1.;

			DVectord v_surf, v_vol_back; 
			T->apply(v_vol, v_surf);
			T->apply_transpose(local_values(local_size(v_surf).get(0), 1.), v_vol_back);


			double min_master = min(v_vol);
			double min_slave = min(v_surf);
			double max_slave = max(v_surf);

			std::cout << "[" << min_slave << ", " << max_slave << "] subset of [" << min_master << ", " << max_master << "]" << std::endl;

			convert(v_surf, *surf_sys.solution);
			surf_sys.solution->close();

			libMesh::Nemesis_IO surf_IO(*surf_mesh);
			surf_IO.write_equation_systems("vol2surf_surf.e", *surf_equation_systems);

			convert(v_vol, *vol_sys.solution);
			convert(v_vol_back, *aux_sys.solution);

			vol_sys.solution->close();
			aux_sys.solution->close();
			libMesh::Nemesis_IO vol_IO(*vol_mesh);
			vol_IO.write_equation_systems("vol2surf_vol.e", *vol_equation_systems);

		} else {
			assert(false);
		}
	}
}