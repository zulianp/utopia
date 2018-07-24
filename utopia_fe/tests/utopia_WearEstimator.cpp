#include "utopia_MechTest.hpp"

#include "utopia_FormEvalTest.hpp"
#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_Mechanics.hpp"
#include "utopia_AffineTransform.hpp"
#include "utopia_Contact.hpp"
#include "utopia_LinearElasticity.hpp"

#include "moonolith_communicator.hpp"
#include "moonolith_describe.hpp"
#include "moonolith_synched_describable.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_refinement.h"


#include "utopia_Socket.hpp"
#include "utopia_NormalTangentialCoordinateSystem.cpp"
#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_ContactSystem.hpp"
#include "utopia_GaitCycle.hpp"
#include "utopia_Wear.hpp"

#include <algorithm>
#include <memory>
#include <array>
#include <fstream>

typedef std::array<double, 2> Point2d;
typedef std::array<double, 3> Point3d;
typedef std::function<std::array<double, 2>(const std::array<double, 2> &p)> Fun2d;
typedef std::function<std::array<double, 3>(const std::array<double, 3> &p)> Fun3d;

namespace utopia {

	class WearEstimator {
	public:

		static void update_wear(
			const double dt,
			const double wear_coefficient,
			const DVectord &sliding_distance,
			const DVectord &normal_stress,
			DVectord &wear
			)
		{	
			wear += (dt * wear_coefficient) * abs(e_mul(sliding_distance, normal_stress));
		}

		void compute_wear_induced_displacement(
			ProductFunctionSpace<LibMeshFunctionSpace> &V,
			const std::vector<int> &boundary_tags,
			DVectord &wear_induced_displacement
			)
		{
			libMesh::MeshBase &mesh = V[0].mesh();
			libMesh::DofMap &dof_map = V[0].dof_map();
			auto dim = mesh.mesh_dimension();

			//compute surface normals
			DVectord is_normal_component;
			DVectord normals;
			DSMatrixd trafo;
			assemble_normal_tangential_transformation(mesh, dof_map, boundary_tags, is_normal_component, normals, trafo);

			//compute surface displacement
			wear_induced_displacement = local_zeros(local_size(wear));
			{
				auto r = range(wear);
				Write<DVectord> w_w(wear_induced_displacement);
				Read<DVectord> r_w(wear), r_n(normals), r_i(is_normal_component);

				for(auto i = r.begin(); i < r.end(); i += dim) {
					if(is_normal_component.get(i) > 0) {
						for(unsigned int d = 0; d < dim; ++d) {
							wear_induced_displacement.set(i + d, normals.get(i + d) * (-extrapolation_factor) * wear.get(i));
						}
					}
				}
			}

		}

		void modify_geometry_with_wear(
			ProductFunctionSpace<LibMeshFunctionSpace> &V,
			const std::vector<int> &boundary_tags)
		{
			libMesh::MeshBase &mesh = V[0].mesh();
			libMesh::DofMap &dof_map = V[0].dof_map();
			auto dim = mesh.mesh_dimension();

			DVectord wear_induced_displacement;
			compute_wear_induced_displacement(V, boundary_tags, wear_induced_displacement);


			// the interior using linear elasticity or laplacian
			// use dirichlet conditions on the whole boundary and warp
			auto u = trial(V);
			auto v = test(V);

			DSMatrixd lapl_mat;
			auto lapl = inner(grad(u), grad(v)) * dX;
			assemble(lapl, lapl_mat);			
			DVectord warped_displacement = local_zeros(local_size(wear_induced_displacement));

			//FIXME warped_displacement passed as dummy
			set_identity_at_constraint_rows(dof_map, lapl_mat);

			Factorization<DSMatrixd, DVectord> solver;
			solver.solve(lapl_mat, wear_induced_displacement, warped_displacement);

			//displace mesh
			apply_displacement(warped_displacement, dof_map, mesh);
			convert(warped_displacement, *V[0].equation_system().solution);
			// convert(wear_induced_displacement, *V[0].equation_system().solution);
			V[0].equation_system().solution->close();

			double wear_magnitude = norm2(wear_induced_displacement);
			std::cout << "wear_magnitude: " << wear_magnitude << std::endl;

			double param_magn = norm2(warped_displacement);
			std::cout << "param_magn: " << param_magn << std::endl;
		}

		void update_aux_system(
			const int main_system_number,
			const MechanicsContext &mech_ctx,
			const MechanicsState &state,
			const Contact &contact,
			libMesh::EquationSystems &es)
		{
			using libMesh::MeshBase;
			auto &main = es.get_system(main_system_number);
			auto &aux = es.get_system("aux");
			const auto &mesh = es.get_mesh();
			const int dim = mesh.mesh_dimension();

			DVectord normal_stress = local_zeros(local_size(state.displacement));
			DVectord sliding_distance = local_zeros(local_size(state.displacement)); 

			// DVectord stress = local_zeros(local_size(state.displacement));

			if(contact.initialized) {
				// stress = e_mul(mech_ctx.inverse_mass_vector, (state.external_force - state.internal_force));
				// stress = (state.external_force - state.internal_force);
				normal_stress = contact.orthogonal_trafo * state.stress;

				DVectord tangential_velocity = contact.orthogonal_trafo * state.velocity;
				sliding_distance = local_zeros(local_size(state.displacement));

				{	
					Read<DVectord>   r_v(tangential_velocity);
					Write<DVectord> w_s(sliding_distance);

					Range r = range(tangential_velocity);

					for(auto i = r.begin(); i < r.end(); i += dim) {

						double dist = 0.;
						for(uint d = 1; d < dim; ++d) {
							auto val = tangential_velocity.get(i + d);
							dist += val*val;
						}

						dist = std::sqrt(dist);
						sliding_distance.set(i, dist);
					}
				}

				sliding_distance = e_mul(contact.is_contact_node, sliding_distance);
				// normal_stress = e_mul(contact.is_contact_node, contact.orthogonal_trafo * stress);

				//override normal stress
				{
					normal_stress *= 0.;

					Read<DVectord> r_s(state.stress), r_n(contact.normals);
					Write<DVectord> w_n(normal_stress);

					Range r = range(state.stress);

					for(auto i = r.begin(); i != r.end(); i+= dim) {
						double ns = 0.;
						for(unsigned int d = 0; d < dim; ++d) {
							ns += state.stress.get(i+d) * contact.normals.get(i+d);
						}

						normal_stress.set(i, ns);
					}
				}

				update_wear(gait_cycle.dt, wear_coefficient, sliding_distance, normal_stress, wear);
				total_wear.push_back(sum(wear));
			}

			{
				Read<DVectord> r_d(state.displacement), r_v(state.velocity);
				Read<DVectord> r_f(state.internal_force), r_ef(state.external_force), r_ns(normal_stress);
				Read<DVectord> r_c(contact.is_contact_node);
				
				auto nd 	= mesh.local_nodes_begin();
				auto nd_end = mesh.local_nodes_end();

				for (; nd != nd_end; ++nd) {
					const libMesh::Node * node = *nd;

					for (unsigned int d = 0; d < dim; ++d) {
						unsigned int source_dof = node->dof_number(main_system_number, d, 0);

						auto dest_dof_disp 	= node->dof_number(aux.number(), var_num_aux[d], 0);
						auto dest_dof_vel 	= node->dof_number(aux.number(), var_num_aux[dim + d], 0);
						auto dest_dof_force = node->dof_number(aux.number(), var_num_aux[2 * dim + d], 0);
						auto dest_dof_ext_force = node->dof_number(aux.number(), var_num_aux[3 * dim + d], 0);

						aux.solution->set(dest_dof_disp,  state.displacement_increment.get(source_dof));
						aux.solution->set(dest_dof_vel,   state.velocity.get(source_dof));
						// aux.solution->set(dest_dof_force, state.internal_force.get(source_dof));

						if(contact.initialized) {
							aux.solution->set(dest_dof_force, state.stress.get(source_dof));
						}

						aux.solution->set(dest_dof_ext_force, state.external_force.get(source_dof));

						auto dest_dof_n = node->dof_number(aux.number(), var_num_aux[4 * dim + d], 0);

						if(contact.initialized) {
							aux.solution->set(dest_dof_n, contact.normals.get(source_dof));
						} else {
							aux.solution->set(dest_dof_n, 0.);
						}

						auto dest_dof_wear_disp = node->dof_number(aux.number(), var_num_aux[5 * dim + d], 0);
						aux.solution->set(dest_dof_wear_disp, wear_induced_displacement.get(source_dof));
					}

					unsigned int source_dof = node->dof_number(main_system_number, 0, 0);
					unsigned int dest_dof_normal_stress = node->dof_number(aux.number(), var_num_aux[6 * dim], 0);
					aux.solution->set(dest_dof_normal_stress, normal_stress.get(source_dof));

					unsigned int dest_dof_sliding_dist = node->dof_number(aux.number(), var_num_aux[6 * dim + 1], 0);
					aux.solution->set(dest_dof_sliding_dist, sliding_distance.get(source_dof));

					unsigned int dest_dof_wear = node->dof_number(aux.number(), var_num_aux[6 * dim + 2], 0);
					aux.solution->set(dest_dof_wear, wear.get(source_dof));

					unsigned int dest_dof_is_contact = node->dof_number(aux.number(), var_num_aux[6 * dim + 3], 0);
					unsigned int dest_dof_is_gap = node->dof_number(aux.number(), var_num_aux[6 * dim + 4], 0);

					if(contact.initialized) {
						aux.solution->set(dest_dof_is_contact, contact.is_contact_node.get(source_dof));
						aux.solution->set(dest_dof_is_gap, contact.gap.get(source_dof));
					} else {
						aux.solution->set(dest_dof_is_contact, 0.);
						aux.solution->set(dest_dof_is_gap, 0.);
					}
				}
			}

			aux.solution->close();
		}

		void init_volume_param_system(libMesh::EquationSystems &es)
		{
			auto &p_sys = es.add_system<libMesh::LinearImplicitSystem>("param");
			param_sys_number = p_sys.number();
		}


		void init_aux_system(
			libMesh::EquationSystems &es,
			libMesh::Order order)
		{
			//init aux system for plotting
			auto &aux = es.add_system<libMesh::LinearImplicitSystem>("aux");
			
			var_num_aux.push_back( aux.add_variable("inc_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("inc_y", libMesh::Order(order), libMesh::LAGRANGE) );

			const int dim = es.get_mesh().mesh_dimension();

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("inc_z", libMesh::Order(order), libMesh::LAGRANGE) );

			var_num_aux.push_back( aux.add_variable("vel_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("vel_y", libMesh::Order(order), libMesh::LAGRANGE) );

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("vel_z", libMesh::Order(order), libMesh::LAGRANGE) );

			var_num_aux.push_back( aux.add_variable("f_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("f_y", libMesh::Order(order), libMesh::LAGRANGE) );

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("f_z", libMesh::Order(order), libMesh::LAGRANGE) );

			var_num_aux.push_back( aux.add_variable("fext_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("fext_y", libMesh::Order(order), libMesh::LAGRANGE) );

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("fext_z", libMesh::Order(order), libMesh::LAGRANGE) );


			var_num_aux.push_back( aux.add_variable("n_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("n_y", libMesh::Order(order), libMesh::LAGRANGE) );

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("n_z", libMesh::Order(order), libMesh::LAGRANGE) );

			var_num_aux.push_back( aux.add_variable("wear_disp_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("wear_disp_y", libMesh::Order(order), libMesh::LAGRANGE) );

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("wear_disp_z", libMesh::Order(order), libMesh::LAGRANGE) );


			var_num_aux.push_back( aux.add_variable("normalstress", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("slidingdistance", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("wear", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("is_contact", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("gap", libMesh::Order(order), libMesh::LAGRANGE) );

			aux.init();	
			aux.update();
		}

		void run(const std::shared_ptr<libMesh::MeshBase> &mesh)
		{
			moonolith::Communicator comm(mesh->comm().get());
			moonolith::root_describe("creating systems....", comm, std::cout);


			const bool override_each_time_step = true;
			const bool is_3d = mesh->mesh_dimension() == 3;

			auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
			auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("wear_test");

			auto elem_order = libMesh::FIRST;
			// auto elem_order = libMesh::SECOND;

			auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
			auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
			auto V  = Vx * Vy;

			if(is_3d) {
				V *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_z");
			}

			auto u  = trial(V);
			auto ux = trial(Vx);
			auto uy = trial(Vy);
			auto vx = test(Vx);
			auto vy = test(Vy);

			auto v = test(V);

			const int dim = mesh->mesh_dimension();


			if(is_3d) {
				auto constr = constraints(
					boundary_conditions(u == coeff(gait_cycle.zero3), {3}),
					boundary_conditions(u == coeff(gait_cycle.bc34),  {4})
				);

				init_constraints(constr);

			} else {
				auto constr = constraints(
					boundary_conditions(u == coeff(gait_cycle.zero2), {3, 4})
				);

				init_constraints(constr);
			}

			Vx.initialize();

			init_aux_system(*equation_systems, elem_order);
			
			//begin: init volume parametrization system
			init_volume_param_system(*equation_systems);

			auto Px = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "u_x", param_sys_number);
			auto Py = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "u_y", param_sys_number);
			auto P = Px * Py;
			auto p = trial(P);

			if(is_3d) {
				P *= LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "u_z", param_sys_number);

				auto p_constr = constraints(
					boundary_conditions(p == coeff(gait_cycle.zero3), {3, 4})
				);

				init_constraints(p_constr);
			} else {
				auto p_constr = constraints(
					boundary_conditions(p == coeff(gait_cycle.zero2), {3, 4})
				);

				init_constraints(p_constr);
			}

			Px.initialize();
			Px.equation_system().update();

			Px.equation_system().solution->zero();
			Px.equation_system().solution->close();


			moonolith::root_describe("DONE", comm, std::cout);
			
			//end: init volume parametrization system

			MechanicsContext mech_ctx;
			mech_ctx.init_mass_matrix(V);

			Size gs({Vx.dof_map().n_dofs()});
			Size ls({Vx.dof_map().n_local_dofs()});
			std::vector<MechanicsState> state(gait_cycle.n_time_steps);
			state[0].init(ls, gs);

			wear = local_zeros(ls);
			wear_induced_displacement = local_zeros(ls);

			mech_ctx.dirichlet_selector = local_values(ls.get(0), 1.);
			apply_zero_boundary_conditions(Vx.dof_map(), mech_ctx.dirichlet_selector);
			mech_ctx.dirichlet_selector = local_values(ls.get(0), 1.) - mech_ctx.dirichlet_selector;

			auto elast = std::make_shared<LinearElasticity<decltype(V), DSMatrixd, DVectord>>(V, params);
			elast->assemble_hessian_and_gradient(state[0].displacement, mech_ctx.stiffness_matrix, state[0].internal_force);
			apply_zero_boundary_conditions(Vx.dof_map(), state[0].external_force);

			auto ef = std::make_shared<ConstantExternalForce>();

			if(is_3d) {
				ef->init((inner(coeff(0.), vx) + inner(coeff(0.), vy) + inner(coeff(0.), v[2])) * dX);
			} else {
				ef->init((inner(coeff(0.), vx) + inner(coeff(0.), vy)) * dX);
			}
			ef->eval(state[0].t, state[0].external_force);
			apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[0].external_force);


			std::cout << "n_dofs: " << Vx.dof_map().n_dofs() << std::endl;

			auto integrator = std::make_shared<ImplicitEuler>(dim, Vx.dof_map());
			// auto cg = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
			// cg->verbose(true);
			// integrator->set_linear_solver(cg);
			
			//begin: set-up semi-geometric multigrid
			// auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
			// auto smoother = std::make_shared<GMRES<DSMatrixd, DVectord> >();
			// auto smg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);
			// smg->verbose(true);
			// smg->init(Vx, 3);
			// integrator->set_linear_solver(smg);
			//end: set-up semi-geometric multigrid

			Contact contact;
			ContactParams contact_params;
			// contact_params.contact_pair_tags = {{2, 1}};
			contact_params.contact_pair_tags = {{1, 2}};
			contact_params.search_radius = 1.7;
			contact_params.use_biorthogonal_basis = static_cast<int>(elem_order) <= 1;

			// libMesh::ExodusII_IO io(*mesh);
			libMesh::Nemesis_IO io(*mesh);


			convert(state[0].displacement, *sys.solution);
			sys.solution->close();
			update_aux_system(0, mech_ctx, state[0], contact, *equation_systems);
			io.write_timestep("wear_test.e", *equation_systems, 1, gait_cycle.t);

			
			

			DVectord overriden_displacement;
			for(std::size_t i = 1; i < gait_cycle.n_time_steps; ++i) {
				std::cout << i << "/" << gait_cycle.n_time_steps << " t = " << gait_cycle.t << std::endl;
				
				overriden_displacement = state[i-1].displacement;

				if(override_each_time_step) {
					gait_cycle.override_displacement(*mesh, Vx.dof_map(), 1, 2, overriden_displacement);
					state[i-1].internal_force *= 0.;
				} 

				gait_cycle.set_time_step(i);
				state[i].init(ls, gs);
				
				//displace mesh
				apply_displacement(overriden_displacement, Vx.dof_map(), *mesh);

				//update boundary conditions				
				if(!contact.init(mesh, make_ref(Vx.dof_map()), contact_params)) {
					std::cerr << "[Error] contact failed" << std::endl;
				}

				contact.print_debug_info();

				state[i].t = gait_cycle.t;				
				ef->eval(state[i].t, state[i].external_force);

				Vx.dof_map().create_dof_constraints(*mesh, gait_cycle.t);
				apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[i].external_force);
				
				//boundary conditions for increment
				// state[i].external_force = e_mul(mech_ctx.dirichlet_selector, state[i].external_force)- e_mul(mech_ctx.dirichlet_selector, state[i-1].displacement);
				if(override_each_time_step) {
					DVectord dummy;
					elast->assemble_hessian_and_gradient(state[i-1].displacement, mech_ctx.stiffness_matrix, dummy);
					double sum_mat = norm2(mech_ctx.stiffness_matrix);
					moonolith::root_describe("stiff_mat_sum: " + std::to_string(sum_mat), comm, std::cout);
				}

				apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[i].external_force);
				state[i].displacement_increment =  e_mul(mech_ctx.dirichlet_selector, state[i].external_force);

				std::cout << "norm(F_e) = " << double(norm2(state[i].external_force)) << std::endl;

				integrator->apply(gait_cycle.dt, mech_ctx, contact, Friction(), state[i-1], state[i]);
				// integrator->apply(gait_cycle.dt, mech_ctx, state[i-1], state[i]);

				// plot_mesh(*mesh, "rot/t_" + std::to_string(i));
				apply_displacement(-overriden_displacement, Vx.dof_map(), *mesh);

				if(override_each_time_step) {
					state[i].displacement = overriden_displacement + state[i].displacement_increment;
					state[i].internal_force = mech_ctx.stiffness_matrix * state[i].displacement_increment;
					apply_zero_boundary_conditions(Vx.dof_map(), state[i].internal_force);
				}

				// state[i].velocity = (1./gait_cycle.dt) * (state[i].displacement - state[i-1].displacement);
				state[i].velocity = (1./gait_cycle.dt) * (state[i].displacement_increment);

				convert(state[i].displacement, *sys.solution);
				sys.solution->close();


				compute_wear_induced_displacement(P, {2}, wear_induced_displacement);
				update_aux_system(0, mech_ctx, state[i], contact, *equation_systems);
				io.write_timestep("wear_test.e", *equation_systems, i + 1, gait_cycle.t);
			}

			// modify_geometry_with_wear(P, {1});
			modify_geometry_with_wear(P, {2});
			io.write_timestep("wear_test.e", *equation_systems, gait_cycle.n_time_steps, gait_cycle.t);

			write_stats();
		}


		void write_stats()
		{
			std::ofstream os("wear_stats.csv");

			if(!os.good()) {
				std::cerr << "[Error] unable to open wear_stats.csv" << std::endl;
				return;
			}


			os << "wear\n";

			for(auto w : total_wear) {
				os << w << "\n";
			}


			os.close();
		}


		WearEstimator()
		: wear_coefficient(7e-3), extrapolation_factor(10.)
		{}

		GaitCycle gait_cycle;
		LameeParameters params;
		std::vector<int> var_num_aux;
		DVectord wear;
		DVectord wear_induced_displacement;
		double wear_coefficient;
		double extrapolation_factor;
		unsigned int param_sys_number;

		std::vector<double> total_wear;
	};

	void run_wear_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		

		// mesh->read("/Users/zulianp/Desktop/algo4u/wearsim/exodus/wear_geoms.e");

		// mesh->read("/Users/zulianp/Desktop/algo4u/wearsim/exodus/knee.e");
		// mesh->read("/Users/zulianp/Desktop/algo4u/wearsim/exodus/knee_fine.e");

		moonolith::Communicator comm(init.comm().get());
		moonolith::root_describe("reading mesh...", comm, std::cout);
		mesh->read("../data/wear_2_far.e");// mesh->all_second_order(false);
		// mesh->read("/Users/zulianp/Desktop/algo4u/wearsim/exodus/toy_coarse.e"); mesh->all_second_order(true);
		moonolith::root_describe("DONE", comm, std::cout);
		// mesh->read("/Users/zulianp/Desktop/algo4u/wearsim/exodus/toy_fine.e");
		

		//libMesh::ExodusII_IO(*mesh).write("exported.e");
		// mesh->read("../data/wear_tri_2.e");

		// unsigned int n_refine = 1;
		// libMesh::MeshRefinement mesh_refinement(*mesh);
		// mesh_refinement.make_flags_parallel_consistent();
		// mesh_refinement.uniformly_refine(n_refine);
		
		WearEstimator we;
		
		//Wear estimator parameters
		we.params = LameeParameters(10., 10.);
		we.params.set_mu(2, 10.);
		we.params.set_lambda(2, 10.);
		we.extrapolation_factor = 1.;

		//gait cycle parameters
		we.gait_cycle.n_time_steps = 100;
		// we.gait_cycle.angle_degree = 45;

		we.gait_cycle.init();


		moonolith::root_describe("running wear simulation....", comm, std::cout);
		we.run(mesh);

		moonolith::root_describe("DONE", comm, std::cout);
	}
}
