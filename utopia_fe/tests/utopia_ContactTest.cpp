#include "utopia_ContactTest.hpp"
#include "utopia_assemble_contact.hpp"


#include "utopia.hpp"

//fe extension
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "ParMortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>

#include "utopia_ContactSimParams.hpp"
#include "libmesh/linear_partitioner.h"
#include "LibmeshContactForMoose.hpp"
#include "LibmeshTransferForMoose.hpp"
#include "LibmeshTransferForMooseReverse.hpp"

#include "utopia_Polygon.hpp"
#include "utopia_NormalTangentialCoordinateSystem.hpp"

#include "moonolith_profiler.hpp"

#include <iostream>

using namespace libMesh;
using std::make_shared;
using std::shared_ptr;

namespace utopia {

	static void solve_contact_problem_2d(
		const libMesh::LibMeshInit &init,
		const std::shared_ptr<libMesh::Mesh> &master_slave,
		const std::vector< std::pair<int, int> > &tags,
		const Real &search_radius
		)
	{
		Chrono c;
		c.start();

		auto order_elem = FIRST;
		int order_quad = order_elem + order_elem;
		int dim = master_slave->mesh_dimension();

		LibMeshFEContext<LinearImplicitSystem> master_slave_context(master_slave);
		auto space_x = fe_space(LAGRANGE, order_elem, master_slave_context);
		auto space_y = fe_space(LAGRANGE, order_elem, master_slave_context);

		auto ux = fe_function(space_x);
		auto uy = fe_function(space_y);
		auto u = prod(ux, uy);

		strong_enforce( boundary_conditions(ux == coeff(0.), {2}) );
		strong_enforce( boundary_conditions(ux == coeff(0.), {4}) );

		strong_enforce( boundary_conditions(uy == coeff(-0.2), {4}) );
		strong_enforce( boundary_conditions(uy == coeff(0.2), {2}) );

		master_slave_context.equation_systems.init();

		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));

		DSMatrixd B;
		moonolith::Communicator moonolith_comm(init.comm().get());

		utopia::DSMatrixd orthogonal_trafos;
		DVectord normals_vec;
		utopia::DVectord gap;
		utopia::DVectord is_contact_node;

		unsigned int variable_number = 0;

		assemble_contact(moonolith_comm, (master_slave), 
			utopia::make_ref(master_slave_context.system.get_dof_map()), 
			variable_number, 
			B, 
			orthogonal_trafos, 
			gap, 
			normals_vec,
			is_contact_node, 
			search_radius,
			tags,
			true);

		DVectord v = local_zeros(local_size(B).get(1));

		each_write(v, [](const SizeType i) -> double {
			return 0.1;
		});

		DVectord mv = B * v;
		DVectord d = sum(B, 1);
		DVectord d_inv = local_zeros(local_size(d));

		{
			Write<DVectord> w_(d_inv);

			each_read(d, [&d_inv](const SizeType i, const double value) {
				if(value < -1e-8) {
					std::cerr << "negative el for " << i << std::endl;
				}

				if(std::abs(value) > 1e-15) {
					d_inv.set(i, 1./value);
				} else {
					d_inv.set(i, 1.);
				}
			});
		}

		DSMatrixd D_inv = diag(d_inv);
		DSMatrixd T = D_inv * B;
		DVectord sum_T = sum(T, 1);

		DVectord D_inv_gap = D_inv * gap;
		T += local_identity(local_size(d).get(0), local_size(d).get(0));

		if(moonolith_comm.is_alone()) plot_scaled_normal_field(*master_slave_context.mesh, normals_vec, D_inv_gap);

		DVectord contact_stress;
		{
			auto b_form = integral(dot(grad(u), grad(u)));

			DenseVector<Real> vec(dim);
			vec.zero();

			auto f 	    = vec_coeff(vec);
			auto l_form = integral(dot(f, u));

			auto ass = make_assembly([&]() -> void {
				double t = MPI_Wtime();

				assemble(u, u, b_form, l_form, *master_slave_context.system.matrix, *master_slave_context.system.rhs);

				t = MPI_Wtime() - t;

				printf("--------------------------------\n");
				printf("Assembly: %g seconds\n", t);
				printf("--------------------------------\n");
			});

			master_slave_context.system.attach_assemble_object(ass);
			master_slave_context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
			master_slave_context.equation_systems.solve();

			DVectord rhs;
			DSMatrixd K;

			convert(*master_slave_context.system.rhs, rhs);
			convert(*master_slave_context.system.matrix, K);

			DVectord sol_c = local_zeros(local_size(rhs));
			DVectord rhs_c = transpose(orthogonal_trafos) * transpose(T) * rhs;
			DSMatrixd K_c  = transpose(orthogonal_trafos) * DSMatrixd(transpose(T) * K * T) * orthogonal_trafos;

			SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<Factorization<DSMatrixd, DVectord> >());
			newton.verbose(true);
			newton.max_it(40);

			newton.set_box_constraints(make_upper_bound_constraints(make_ref(D_inv_gap)));
			newton.solve(K_c, rhs_c, sol_c);

			DVectord sol = T * (orthogonal_trafos * sol_c);
			convert(sol, *master_slave_context.system.solution);
		}

		ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("sol2.e", master_slave_context.equation_systems);

		convert(is_contact_node, *master_slave_context.system.solution);
		ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("is_c_node2.e", master_slave_context.equation_systems);

				// convert(gap, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("gap.e", master_slave_context.equation_systems);

				// convert(normals_vec, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("normals.e", master_slave_context.equation_systems);

				// convert(d, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("d.e", master_slave_context.equation_systems);

				// normals_vec = orthogonal_trafos * normals_vec;
				// convert(normals_vec, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("H_n.e", master_slave_context.equation_systems);
	}


	static void solve_contact_problem_3d(
		const libMesh::LibMeshInit &init,
		const std::shared_ptr<libMesh::Mesh> &master_slave,
		const std::vector< std::pair<int, int> > &tags,
		const Real &search_radius
		)
	{
		Chrono c;
		c.start();

		auto order_elem = FIRST;
		int order_quad = order_elem + order_elem;
		int dim = master_slave->mesh_dimension();

		LibMeshFEContext<LinearImplicitSystem> master_slave_context(master_slave);
		auto space_x = fe_space(LAGRANGE, order_elem, master_slave_context);
		auto space_y = fe_space(LAGRANGE, order_elem, master_slave_context);
		auto space_z = fe_space(LAGRANGE, order_elem, master_slave_context);

		auto ux = fe_function(space_x);
		auto uy = fe_function(space_y);
		auto uz = fe_function(space_z);
		auto u = prod(ux, uy, uz);

		strong_enforce( boundary_conditions(ux == coeff(0.), {2}) );
		strong_enforce( boundary_conditions(ux == coeff(0.), {4}) );

		strong_enforce( boundary_conditions(uy == coeff(0.), {4}) );
		strong_enforce( boundary_conditions(uy == coeff(0.), {2}) );

		strong_enforce( boundary_conditions(uz == coeff(1.), {4}) );
		strong_enforce( boundary_conditions(uz == coeff(-1.), {2}) );

		master_slave_context.equation_systems.init();

		ux.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uy.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));
		uz.set_quad_rule(make_shared<libMesh::QGauss>(dim, SECOND));

		DSMatrixd B;
		moonolith::Communicator moonolith_comm(init.comm().get());

		utopia::DSMatrixd orthogonal_trafos;
		DVectord normals_vec;
		utopia::DVectord gap;
		utopia::DVectord is_contact_node;

		unsigned int variable_number = 0;

		assemble_contact(moonolith_comm, (master_slave), 
			utopia::make_ref(master_slave_context.system.get_dof_map()), 
			variable_number, 
			B, 
			orthogonal_trafos, 
			gap, 
			normals_vec,
			is_contact_node, 
			search_radius,
			tags,
			true);

		DVectord v = local_zeros(local_size(B).get(1));

		each_write(v, [](const SizeType i) -> double {
			return 0.1;
		});

		DVectord mv = B * v;
		DVectord d = sum(B, 1);
		DVectord d_inv = local_zeros(local_size(d));

		{
			Write<DVectord> w_(d_inv);

			each_read(d, [&d_inv](const SizeType i, const double value) {
				if(value < -1e-8) {
					std::cerr << "negative el for " << i << std::endl;
				}

				if(std::abs(value) > 1e-15) {
					d_inv.set(i, 1./value);
				} else {
					d_inv.set(i, 1.);
				}
			});
		}

		DSMatrixd D_inv = diag(d_inv);
		DSMatrixd T = D_inv * B;
		DVectord sum_T = sum(T, 1);

		DVectord D_inv_gap = D_inv * gap;
		T += local_identity(local_size(d).get(0), local_size(d).get(0));

		if(moonolith_comm.is_alone()) plot_scaled_normal_field(*master_slave_context.mesh, normals_vec, D_inv_gap);

		DVectord contact_stress;
		{
			auto b_form = integral(dot(grad(u), grad(u)));

			DenseVector<Real> vec(dim);
			vec.zero();

			auto f 	    = vec_coeff(vec);
			auto l_form = integral(dot(f, u));

			auto ass = make_assembly([&]() -> void {
				double t = MPI_Wtime();

				assemble(u, u, b_form, l_form, *master_slave_context.system.matrix, *master_slave_context.system.rhs);

				t = MPI_Wtime() - t;

				printf("--------------------------------\n");
				printf("Assembly: %g seconds\n", t);
				printf("--------------------------------\n");
			});

			master_slave_context.system.attach_assemble_object(ass);
			master_slave_context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
			master_slave_context.equation_systems.solve();

			DVectord rhs;
			DSMatrixd K;

			convert(*master_slave_context.system.rhs, rhs);
			convert(*master_slave_context.system.matrix, K);

			DVectord sol_c = local_zeros(local_size(rhs));
			DVectord rhs_c = transpose(orthogonal_trafos) * transpose(T) * rhs;
			DSMatrixd K_c  = transpose(orthogonal_trafos) * DSMatrixd(transpose(T) * K * T) * orthogonal_trafos;

			SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<Factorization<DSMatrixd, DVectord> >());
			newton.verbose(true);
			newton.max_it(40);

			newton.set_box_constraints(make_upper_bound_constraints(make_ref(D_inv_gap)));
			newton.solve(K_c, rhs_c, sol_c);

			DVectord sol = T * (orthogonal_trafos * sol_c);
			convert(sol, *master_slave_context.system.solution);
		}

		ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("sol3.e", master_slave_context.equation_systems);

		convert(is_contact_node, *master_slave_context.system.solution);
		ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("is_c_node3.e", master_slave_context.equation_systems);

				// convert(gap, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("gap.e", master_slave_context.equation_systems);

				// convert(normals_vec, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("normals.e", master_slave_context.equation_systems);

				// convert(d, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("d.e", master_slave_context.equation_systems);

				// normals_vec = orthogonal_trafos * normals_vec;
				// convert(normals_vec, *master_slave_context.system.solution);
				// ExodusII_IO(*master_slave_context.mesh).write_equation_systems ("H_n.e", master_slave_context.equation_systems);
	}


	void run_contact_test(LibMeshInit &init)
	{
		// auto mesh = make_shared<Mesh>(init.comm());
		// mesh->read("../data/fine_contact_2d.e");
		// // mesh->read("../data/hertz_2d.e");
		// Real search_radius = 0.1;
		// solve_contact_problem_2d(init, mesh, {{102, 101}}, search_radius);


		auto mesh = make_shared<Mesh>(init.comm());
		// mesh->read("../data/hertz_530.e");
		// mesh->read("../data/quasi_signorini_4593.e");
		mesh->read("../data/two_rocks_26653.e");
		Real search_radius = 1.;
		solve_contact_problem_3d(init, mesh, {{1, 3}}, search_radius);
	}
}