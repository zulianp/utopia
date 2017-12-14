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

#include "moonolith_communicator.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/reference_counted_object.h"

#include "utopia_Socket.hpp"
#include "utopia_NormalTangentialCoordinateSystem.cpp"

#include <algorithm>
#include <memory>
#include <array>

typedef std::array<double, 2> Point2d;
typedef std::array<double, 3> Point3d;
typedef std::function<std::array<double, 2>(const std::array<double, 2> &p)> Fun2d;

namespace utopia {

	// static bool assemble_normals_on_linear_surface(
	// 	const libMesh::MeshBase &mesh,
	// 	const libMesh::DofMap &dof_map,
	// 	const std::vector<int> &boundary_tags,
	// 	DVectord &is_normal_component,
	// 	DVectord &normals)
	// {
	//     const SizeType n_local_dofs = dof_map.n_local_dofs();
	//     const unsigned int dim = mesh.mesh_dimension();

	//     DVectord normals_buff = local_zeros(n_local_dofs);


	//     QGauss quad(dim-1, SECOND);
	//     std::unique_ptr<FEBase> fe = FEBase::build(n_dims, FIRST);
	//     const auto & phi = fe->get_phi();
	//     const auto & JxW = fe->get_JxW();

	//     Point normal;
	//     std::vector<dof_id_type> dof_indices;

	//     for(auto e_it = mesh.active_local_elements_begin();
	//         e_it != mesh.active_local_elements_end(); ++e_it)
	//     {
	//         const auto &e = **e_it;
	//         dof_map.dof_indices(&e, dof_indices);

	//         for(uint side = 0; side < e.n_sides(); ++side) {
	//             if(e.neighbor_ptr(side) != nullptr) { continue; }
	//              compute_side_normal(dim, *side, normal);


	//              for(uint i = 0; i < n_test; ++i) {
	//              	for(uint qp = 0; qp < n_qp; ++qp) {

	//              		for(uint d = 0; d < dim; ++d) {
	//              			normal_vec(i, d) += phi[i][qp] * normal(d) * JxW[qp];
	//              		}
	//              	}
	//              }

	//         }
	//     }
	// }



	static void apply_displacement(
		const DVectord &displacement_increment,
		const libMesh::DofMap &dof_map,
		libMesh::MeshBase &mesh)
	{	
		//FIXME
		int sys_num = 0;

		moonolith::Communicator comm(mesh.comm().get());

		if(!comm.is_alone()) {
			auto r = range(displacement_increment);
			Read<DVectord> r_d(displacement_increment);

			auto m_begin = mesh.active_local_elements_begin();
			auto m_end   = mesh.active_local_elements_end();

			std::vector<PetscInt> idx;
			std::set<PetscInt> unique_idx;
			std::map<libMesh::dof_id_type, double> idx_to_value;
			std::vector<libMesh::dof_id_type> dof_indices;
			
			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				dof_map.dof_indices(*m_it, dof_indices);
				for(auto dof_id : dof_indices) {
					if(r.inside(dof_id)) {
						idx_to_value[dof_id] = displacement_increment.get(dof_id);
					} else {
						unique_idx.insert(dof_id);
					}
				}
			}

			idx.insert(idx.end(), unique_idx.begin(), unique_idx.end());
			DVectord out = displacement_increment.select(idx);
			{
				Read<DVectord> r_out(out);
				auto range_out = range(out);

				for(std::size_t i = 0; i < idx.size(); ++i) {
					idx_to_value[idx[i]] = out.get(range_out.begin() + i);
				}
			}	

			for(auto m_it = m_begin; m_it != m_end; ++m_it) { 
				auto &e = **m_it;
				for(int i = 0; i < e.n_nodes(); ++i) {
					auto &node = e.node_ref(i);

					for(unsigned int c = 0; c < mesh.mesh_dimension(); ++c) {
						const int dof_id = node.dof_number(sys_num, c, 0);
						assert(idx_to_value.find(dof_id) != idx_to_value.end());
						double &val = idx_to_value[dof_id];
						node(c) += val;
						val = 0.;
					}
				}
			}

		} else {

			Read<DVectord> r_d(displacement_increment);

			auto m_it  = mesh.local_nodes_begin();
			auto m_end = mesh.local_nodes_end();


			for(; m_it != m_end; ++m_it) { 
				for(unsigned int c = 0; c < mesh.mesh_dimension(); ++c) {
					const int dof_id = (*m_it)->dof_number(sys_num, c, 0);
					(**m_it)(c) += displacement_increment.get(dof_id);
				}
			}
		}
	}


	class GaitCycle {
	public:

		GaitCycle()
		{
			n_time_steps = 50;
			t_end = 10.;
			angle_degree = 60.;
			init();
		}

		void set_time_step(const std::size_t time_step)
		{
			t = time_step * dt;
			if(time_step > n_time_steps/2) {
				negative_dir = true;
			}
		}

		void toggle_dir()
		{
			negative_dir = !negative_dir;
		}

		void init()
		{
			t = 0.;
			dt = (t_end - t)/(n_time_steps - 1);	
			angle_radian = (angle_degree/180 * M_PI);
			d_angle = angle_radian/(t_end - t);
			negative_dir = false;

			rotate2 = [this](const Point2d &p) -> Point2d {
				AffineTransform trafo;
				trafo.make_rotation(2, this->t * this->d_angle, 'y');
				trafo.translation[0] = -p[0];
				trafo.translation[1] = -p[1];
				return trafo.apply(p);
			};

			zero2 = [](const Point2d &p) -> Point2d {
				return {0., 0.};
			};

			translate2_y = [this](const Point2d &p) -> Point2d {
				return {0., std::min(2*this->dt*0.1, 0.02) };
			};
		}

		void override_displacement(
			const libMesh::MeshBase &mesh,
			const libMesh::DofMap &dof_map,
			const int block_id_rot,
			const int block_id_trasl,
			DVectord &displacement) const
		{
			//FIXME
			const int main_system_number = 0;
			// displacement = local_zeros(dof_map.n_local_dofs());

			const auto dim = mesh.mesh_dimension();

			auto r = range(displacement);
			Write<DVectord> w_d(displacement);

			for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
				const auto &e = **e_it;
				if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

				for(std::size_t i = 0; i < e.n_nodes(); ++i) {
					const auto &node = e.node_ref(i);

					Point2d p{ node(0), node(1) };

					if(e.subdomain_id() == block_id_rot){
						p = rotate2(p);
					} else if(e.subdomain_id() == block_id_trasl) {
						p = translate2_y(p);
					}

					for(unsigned int d = 0; d < dim; ++d) {
						unsigned int dof = node.dof_number(main_system_number, d, 0);

						if(r.inside(dof)) {
							displacement.set(dof, p[d]);
						}
					}
				}
			}
		}

		std::size_t n_time_steps;
		double dt;
		double t;
		double t_end;
		double angle_degree;
		double angle_radian;
		double d_angle;
		Fun2d rotate2;
		Fun2d translate2_y;
		Fun2d zero2;
		bool negative_dir;
	};

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

		void modify_geometry_with_wear(
			ProductFunctionSpace<LibMeshFunctionSpace> &V,
			const std::vector<int> &boundary_tags)
		{
			libMesh::MeshBase &mesh = V[0].mesh();
			libMesh::DofMap &dof_map = V[0].dof_map();
			auto dim = mesh.mesh_dimension();

			//compute surface normals
			DVectord is_normal_component;
			DVectord normals;
			DSMatrixd trafo;

			//why this does not work????
			assemble_normal_tangential_transformation(mesh, dof_map, boundary_tags, is_normal_component, normals, trafo);

			//compute surface displacement
			DVectord wear_induced_displacement = local_zeros(local_size(wear));
			{
				auto r = range(wear);
				Write<DVectord> w_w(wear_induced_displacement);
				Read<DVectord> r_w(wear), r_n(normals), r_i(is_normal_component);

				for(auto i = r.begin(); i != r.end(); i+= dim) {
					if(is_normal_component.get(i) > 0) {
						for(unsigned int d = 0; d < dim; ++d) {
							wear_induced_displacement.set(i + d, normals.get(i + d) * (-extrapolation_factor) * wear.get(i));
							// wear_induced_displacement.set(i + d, normals.get(i + d) * 0.1);
						}
					}
				}
			}
			
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

			if(contact.initialized) {
				DVectord stress = e_mul(mech_ctx.inverse_mass_vector, (state.external_force - state.internal_force));
				apply_zero_boundary_conditions(main.get_dof_map(), stress);

				normal_stress = contact.orthogonal_trafo * stress;

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
							dist = val*val;
						}

						dist = std::sqrt(dist);
						sliding_distance.set(i, dist);
					}
				}

				update_wear(gait_cycle.dt, wear_coefficient, sliding_distance, normal_stress, wear);
				sliding_distance = e_mul(contact.is_contact_node, sliding_distance);
				normal_stress = e_mul(contact.is_contact_node, contact.orthogonal_trafo * stress);
			}

			{
				Read<DVectord> r_d(state.displacement), r_v(state.velocity);
				Read<DVectord> r_f(state.internal_force), r_ef(state.external_force), r_ns(normal_stress);
				
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
						aux.solution->set(dest_dof_force, state.internal_force.get(source_dof));
						aux.solution->set(dest_dof_ext_force, state.external_force.get(source_dof));
					}

					unsigned int source_dof = node->dof_number(main_system_number, 0, 0);
					unsigned int dest_dof_normal_stress = node->dof_number(aux.number(), var_num_aux[4 * dim], 0);
					aux.solution->set(dest_dof_normal_stress, normal_stress.get(source_dof));

					unsigned int dest_dof_sliding_dist = node->dof_number(aux.number(), var_num_aux[4 * dim + 1], 0);
					aux.solution->set(dest_dof_sliding_dist, sliding_distance.get(source_dof));

					unsigned int dest_dof_wear = node->dof_number(aux.number(), var_num_aux[4 * dim + 2], 0);
					aux.solution->set(dest_dof_wear, wear.get(source_dof));
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
			const int order = 1)
		{
			//init aux system for plotting
			auto &aux = es.add_system<libMesh::LinearImplicitSystem>("aux");
			
			var_num_aux.push_back( aux.add_variable("inc_x", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("inc_y", libMesh::Order(order), libMesh::LAGRANGE) );

			const int dim = es.get_mesh().mesh_dimension();

			if(dim > 2) 
				var_num_aux.push_back( aux.add_variable("disp_z", libMesh::Order(order), libMesh::LAGRANGE) );

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

			var_num_aux.push_back( aux.add_variable("normalstress", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("slidingdistance", libMesh::Order(order), libMesh::LAGRANGE) );
			var_num_aux.push_back( aux.add_variable("wear", libMesh::Order(order), libMesh::LAGRANGE) );

			aux.init();	
			aux.update();
		}

		void run(const std::shared_ptr<libMesh::MeshBase> &mesh)
		{
			const bool override_each_time_step = true;

			auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
			auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("wear_test");

			auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "disp_x");
			auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "disp_y");
			auto V  = Vx * Vy;

			auto u  = trial(V);
			auto ux = trial(Vx);
			auto uy = trial(Vy);
			auto vx = test(Vx);
			auto vy = test(Vy);

			const int dim = mesh->mesh_dimension();

			auto constr = constraints(
				// boundary_conditions(u == coeff(gait_cycle.translate2_y), {4}),
				boundary_conditions(u == coeff(gait_cycle.zero2),        {3, 4})
			);

			FEBackend<LIBMESH_TAG>::init_constraints(constr);
			Vx.initialize();

			init_aux_system(*equation_systems, 1);
			
			//begin: init volume parametrization system
			init_volume_param_system(*equation_systems);

			auto Px = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_x", param_sys_number);
			auto Py = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u_y", param_sys_number);
			auto P = Px * Py;

			auto p_constr = constraints(
				boundary_conditions(trial(P) == coeff(gait_cycle.zero2), {1, 2, 3, 4})
			);

			FEBackend<LIBMESH_TAG>::init_constraints(p_constr);
			Px.initialize();
			Px.equation_system().update();

			Px.equation_system().solution->zero();
			Px.equation_system().solution->close();

			//end: init volume parametrization system

			MechanicsContext mech_ctx;
			mech_ctx.init_mass_matrix(V);

			Size gs({Vx.dof_map().n_dofs()});
			Size ls({Vx.dof_map().n_local_dofs()});
			std::vector<MechanicsState> state(gait_cycle.n_time_steps);
			state[0].init(ls, gs);

			wear = local_zeros(ls);

			mech_ctx.dirichlet_selector = local_values(ls.get(0), 1.);
			apply_zero_boundary_conditions(Vx.dof_map(), mech_ctx.dirichlet_selector);
			mech_ctx.dirichlet_selector = local_values(ls.get(0), 1.) - mech_ctx.dirichlet_selector;

			auto elast = std::make_shared<LinearElasticity>();
			elast->init(V, params, state[0].displacement, mech_ctx.stiffness_matrix, state[0].internal_force);
			apply_zero_boundary_conditions(Vx.dof_map(), state[0].external_force);

			auto ef = std::make_shared<ConstantExternalForce>();
			ef->init((inner(coeff(0.), vx) + inner(coeff(0.), vy)) * dX);
			ef->eval(state[0].t, state[0].external_force);
			apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[0].external_force);

			auto integrator = std::make_shared<ImplicitEuler>(dim, Vx.dof_map());

			Contact contact;
			ContactParams contact_params;
			contact_params.contact_pair_tags = {{2, 1}};
			contact_params.search_radius = 1.;

			// libMesh::ExodusII_IO io(*mesh);
			libMesh::Nemesis_IO io(*mesh);

			convert(state[0].displacement, *sys.solution);
			sys.solution->close();
			update_aux_system(0, mech_ctx, state[0], contact, *equation_systems);
			io.write_timestep("wear_test.e", *equation_systems, 1, gait_cycle.t);

			
			DVectord overriden_displacement;
			for(std::size_t i = 1; i < gait_cycle.n_time_steps; ++i) {
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

				state[i].t = gait_cycle.t;				
				ef->eval(state[i].t, state[i].external_force);

				Vx.dof_map().create_dof_constraints(*mesh, gait_cycle.t);
				apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[i].external_force);
				
				//boundary conditions for increment
				// state[i].external_force = e_mul(mech_ctx.dirichlet_selector, state[i].external_force)- e_mul(mech_ctx.dirichlet_selector, state[i-1].displacement);
				if(override_each_time_step) {
					elast->assemble_hessian(V, params, state[i-1].displacement, mech_ctx.stiffness_matrix);
				}

				apply_boundary_conditions(Vx.dof_map(), mech_ctx.stiffness_matrix, state[i].external_force);
				integrator->apply(gait_cycle.dt, mech_ctx, contact, Friction(), state[i-1], state[i]);
				// integrator->apply(gait_cycle.dt, mech_ctx, state[i-1], state[i]);

				// plot_mesh(*mesh, "rot/t_" + std::to_string(i));
				apply_displacement(-overriden_displacement, Vx.dof_map(), *mesh);

				if(override_each_time_step) {
					state[i].displacement = overriden_displacement + state[i].displacement_increment;
					state[i].internal_force = mech_ctx.stiffness_matrix * state[i].displacement_increment;
					apply_zero_boundary_conditions(Vx.dof_map(), state[i].internal_force);
				}

				state[i].velocity = (1./gait_cycle.dt) * (state[i].displacement - state[i-1].displacement);

				convert(state[i].displacement, *sys.solution);
				sys.solution->close();

				update_aux_system(0, mech_ctx, state[i], contact, *equation_systems);
				io.write_timestep("wear_test.e", *equation_systems, i + 1, gait_cycle.t);
			}

			modify_geometry_with_wear(P, {1});
			io.write_timestep("wear_test.e", *equation_systems, gait_cycle.n_time_steps, gait_cycle.t);
			plot_mesh(*mesh, "worn_geometry");

		}


		WearEstimator()
		: wear_coefficient(7e-3), extrapolation_factor(10.)
		{}

		GaitCycle gait_cycle;
		LameeParameters params;
		std::vector<int> var_num_aux;
		DVectord wear;
		double wear_coefficient;
		double extrapolation_factor;
		unsigned int param_sys_number;

	};

	void run_wear_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());		
		mesh->read("../data/wear_2.e");
		// mesh->read("../data/wear_tri_2.e");

		
		WearEstimator we;
		
		//Wear estimator parameters
		we.params = LameeParameters(10., 10.);
		we.params.set_mu(2, 1.);
		we.params.set_lambda(2, 1.);
		we.extrapolation_factor = 10.;

		//gait cycle parameters
		we.gait_cycle.n_time_steps = 100;
		// we.gait_cycle.angle_degree = 10;
		we.gait_cycle.init();

		we.run(mesh);
	}
}
