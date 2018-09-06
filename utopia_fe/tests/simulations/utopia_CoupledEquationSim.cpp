#include "utopia_CoupledEquationSim.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_MeshTransferOperator.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_TransferAssembler.hpp"

#include "libmesh/mesh_refinement.h"
#include <memory>

namespace utopia {

	typedef utopia::LibMeshFunctionSpace FunctionSpaceT;


	bool assemble_projection(FunctionSpaceT &from, FunctionSpaceT &to, DSMatrixd &B, DSMatrixd &D)
	{
		auto assembler = std::make_shared<L2LocalAssembler>(from.mesh().mesh_dimension(), false, true);
		auto local2global = std::make_shared<Local2Global>(false);

		TransferAssembler transfer_assembler(assembler, local2global);

		std::vector< std::shared_ptr<DSMatrixd> > mats;
		if(!transfer_assembler.assemble(
			make_ref(from.mesh()),
			make_ref(from.dof_map()),
			make_ref(to.mesh()),
			make_ref(to.dof_map()),
			mats)) {
			return false;
		}

		B = std::move(*mats[0]);
		D = std::move(*mats[1]);
		return true;
	}

	static void refine(const int n_refs, libMesh::MeshBase &mesh)
	{
		if(n_refs <= 0) return;

		libMesh::MeshRefinement mesh_refinement(mesh);
		mesh_refinement.make_flags_parallel_consistent();
		mesh_refinement.uniformly_refine(n_refs);
	}

	void run_coupled_equation_test(libMesh::LibMeshInit &init)
	{
		

		//model parameters
		const unsigned int n = 100;
		const unsigned int m = 10;

		//discretization parameters
		const auto elem_type = libMesh::QUAD8;
		const auto elem_order = libMesh::FIRST;

		auto mesh_master = std::make_shared<libMesh::DistributedMesh>(init.comm());
		libMesh::MeshTools::Generation::build_square(
			*mesh_master,
			n, n,
			0, 1.,
			0, 1.,
			elem_type
		);

		auto mesh_slave = std::make_shared<libMesh::DistributedMesh>(init.comm());
		// libMesh::MeshTools::Generation::build_square(
		// 	*mesh_slave,
		// 	m, m/2,
		// 	0., 1.,
		// 	0.35, 0.4,
		// 	elem_type
		// );

		mesh_slave->read("../data/frac/line.e");
		// refine(1, *mesh_slave);


		//equations system
		auto equation_systems_master = std::make_shared<libMesh::EquationSystems>(*mesh_master);
		auto &sys_master = equation_systems_master->add_system<libMesh::LinearImplicitSystem>("master");

		auto equation_systems_slave = std::make_shared<libMesh::EquationSystems>(*mesh_slave);
		auto &sys_slave = equation_systems_slave->add_system<libMesh::LinearImplicitSystem>("slave");

		//scalar function space
		auto V_m = FunctionSpaceT(equation_systems_master, libMesh::LAGRANGE, elem_order, "u");
		auto V_s = FunctionSpaceT(equation_systems_slave, libMesh::LAGRANGE, elem_order, "u");


		auto u_m = trial(V_m);
		auto v_m = test(V_m);

		auto u_s = trial(V_s);
		auto v_s = test(V_s);

		init_constraints(
			constraints(
				boundary_conditions(u_m == coeff(4.), {1}),
				boundary_conditions(u_m == coeff(1.), {3})
			)
		);

		init_constraints(
			constraints(
				boundary_conditions(u_s == coeff(1.), {1}),
				boundary_conditions(u_s == coeff(4.), {2})
			)
		);


		V_m.initialize();
		V_s.initialize();

		auto eq_m = inner(grad(u_m), grad(v_m)) * dX == inner(coeff(0.), v_m) * dX;
		auto eq_s = 10. * inner(grad(u_s), grad(v_s)) * dX == inner(coeff(0.), v_m) * dX;

		DSMatrixd A_m, A_s;
		DVectord rhs_m, rhs_s;
		utopia::assemble(eq_m, A_m, rhs_m);
		utopia::assemble(eq_s, A_s, rhs_s);

		apply_boundary_conditions(V_m.dof_map(), A_m, rhs_m);
		apply_boundary_conditions(V_s.dof_map(), A_s, rhs_s);


		// disp(A_m);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		DSMatrixd B, D;
		assemble_projection(V_m, V_s, B, D);

		D *= -1.;
		// B *= -1.;

		auto s_m = local_size(A_m);
		auto s_s = local_size(A_s);


		auto nnz_x_row_m = 
		std::max(*std::max_element(V_m.dof_map().get_n_nz().begin(), V_m.dof_map().get_n_nz().end()),
			     *std::max_element(V_m.dof_map().get_n_oz().begin(), V_m.dof_map().get_n_oz().end()));

		auto nnz_x_row_s = 
		std::max(*std::max_element(V_s.dof_map().get_n_nz().begin(), V_s.dof_map().get_n_nz().end()),
			     *std::max_element(V_s.dof_map().get_n_oz().begin(), V_s.dof_map().get_n_oz().end()));


		DSMatrixd A = local_sparse(
			s_m.get(0) + 2 * s_s.get(0),
			s_m.get(1) + 2 * s_s.get(1),
			2*nnz_x_row_m + 2*nnz_x_row_s //FIXME
		);


		auto rr_m = row_range(A_m);
		auto rr_s = row_range(A_s);

		auto off_r   = size(A_m).get(0);
		auto off_c   = size(A_m).get(1);

		auto off_r_l = off_r + size(A_s).get(0);
		auto off_c_l = off_c + size(A_s).get(1);

		{
			Write<DSMatrixd> w_(A);

			for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
				{
					RowView<DSMatrixd> row(A_m, i);

					for(auto k = 0; k < row.n_values(); ++k) {
						A.set(i, row.col(k), row.get(k));
					}
				}
			}

			for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
				
				{
					RowView<DSMatrixd> row(A_s, i);

					for(auto k = 0; k < row.n_values(); ++k) {
						A.set(off_r + i, off_c + row.col(k), row.get(k));
					}
				}

				{
					RowView<DSMatrixd> row(B, i);

					for(auto k = 0; k < row.n_values(); ++k) {
						//(1, 3) 
						if(!V_m.dof_map().is_constrained_dof(row.col(k))) {
							A.set(row.col(k), off_c_l + i, row.get(k));
						}

						//(3, 1)
						A.set(off_r_l + i, row.col(k), row.get(k));
					}
				}

				{
					RowView<DSMatrixd> row(D, i);

					for(auto k = 0; k < row.n_values(); ++k) {

						//(2, 3)
						if(!V_s.dof_map().is_constrained_dof(row.col(k))) {
							A.set(off_r + row.col(k), off_c_l + i, row.get(k));
						}

						A.set(off_r_l + i, off_c + row.col(k), row.get(k));
					}
				}

			}
		}

		DVectord rhs = local_zeros(local_size(rhs_m).get(0) + 2*local_size(rhs_s).get(0));

		{
			Write<DVectord> w_(rhs);
			Write<DVectord> r_m(rhs_m), r_s(rhs_s);

			for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
				rhs.set(i, rhs_m.get(i));				
			}

			for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
				rhs.set(i, off_r + rhs_s.get(i));
			}
		}



		A.implementation().set_name("a");
		write("A.m", A);

		// disp(A);

		Factorization<DSMatrixd, DVectord> op;
		// GMRES<DSMatrixd, DVectord> op;
		// op.max_it(300000);
		// op.verbose(true);
		op.update(make_ref(A));


		DVectord sol = local_zeros(local_size(rhs));
		op.apply(rhs, sol);


		DVectord sol_m  = local_zeros(local_size(rhs_m));
		DVectord sol_s  = local_zeros(local_size(rhs_s));
		DVectord lagr   = local_zeros(local_size(rhs_s));

		{
			Write<DVectord> w_(sol);
			Write<DVectord> r_m(sol_m), r_s(sol_s);

			for(auto i = rr_m.begin(); i != rr_m.end(); ++i) {
				sol_m.set(i, sol.get(i));				
			}

			for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
				sol_s.set(i, sol.get(i) - off_r);
			}

			for(auto i = rr_s.begin(); i != rr_s.end(); ++i) {
				lagr.set(i, sol.get(i) - off_r_l);
			}
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Factorization<DSMatrixd, DVectord> op_m;
		// op_m.update(make_ref(A_m));


		// Factorization<DSMatrixd, DVectord> op_s;
		// op_s.update(make_ref(A_s));


		// DVectord lagr_m = local_zeros(local_size(rhs_m));
		// DVectord lagr_s = local_zeros(local_size(rhs_s));



		// DVectord rhs_lagr_m;
		// DVectord delta_lagr;

		// double dumping = 1.;


		// MeshTransferOperator t(mesh_master,
		// 					   make_ref(V_m.dof_map()),
		// 					   mesh_slave,
		// 					   make_ref(V_s.dof_map())
		// 					   );


		// // t.initialize(INTERPOLATION);
		// t.initialize(L2_PROJECTION);


		// for(int i = 0; i < 20; ++i) {
		// 	apply_zero_boundary_conditions(V_m.dof_map(), lagr_m);
		// 	rhs_lagr_m = rhs_m + lagr_m;
			
		// 	op_m.apply(rhs_lagr_m, sol_m);

		// 	t.apply(sol_m, sol_s);

		// 	lagr_s = rhs_s - A_s * sol_s;

		// 	double n_lagr_s = norm2(lagr_s);

		// 	disp(n_lagr_s);

		// 	if(n_lagr_s < 1e-8) {
		// 		break;
		// 	}

		// 	t.apply_transpose(lagr_s, delta_lagr);
		// 	lagr_m += dumping * delta_lagr;
		// }


		// sol_m.set(0.);
		// sol_s.set(0.);
		// apply_boundary_conditions(V_m.dof_map(), A_m, sol_m);
		// apply_boundary_conditions(V_s.dof_map(), A_s, sol_s);


		libMesh::ExodusII_IO io_m(*mesh_master);
		libMesh::ExodusII_IO io_s(*mesh_slave);
		libMesh::ExodusII_IO io_l(*mesh_slave);

		utopia::convert(sol_m, *V_m.equation_system().solution);
		V_m.equation_system().solution->close();
		io_m.write_timestep(V_m.equation_system().name() + ".e", V_m.equation_systems(), 1, 0);


		utopia::convert(sol_s, *V_s.equation_system().solution);
		V_s.equation_system().solution->close();
		io_s.write_timestep(V_s.equation_system().name() + ".e", V_s.equation_systems(), 1, 0);

		utopia::convert(lagr, *V_s.equation_system().solution);
		V_s.equation_system().solution->close();
		io_l.write_timestep("lagr.e", V_s.equation_systems(), 1, 0);
	}
}
