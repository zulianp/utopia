#include "utopia_VolumeInterpolationTest.hpp"

#include "utopia_libmesh.hpp"
#include "moonolith_communicator.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_TransferAssembler.hpp"

#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/mesh_refinement.h"
#include <algorithm>
#include <cmath>

typedef utopia::LibMeshFunctionSpace FunctionSpaceT;

namespace utopia {

	void run_volume_interpolation_test(libMesh::LibMeshInit &init)
	{
		auto n_master = 4;
		auto n_slave  = 5;
		auto elem_order_master = libMesh::SECOND;
		auto elem_order_slave  = libMesh::SECOND;

		// auto elem_order_master = libMesh::FIRST;
		// auto elem_order_slave  = libMesh::FIRST;

		auto master_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());
		auto slave_mesh  = std::make_shared<libMesh::DistributedMesh>(init.comm());

		bool test_files = true;
		bool two_d      = false;

		if(test_files) {
			std::string path = "../data/test/quad_1.e";
			// std::string path = "../data/test/LV_ellipsoid_tet_Z.e";

			master_mesh->read(path);
			slave_mesh->read(path);

			master_mesh->all_second_order(false);
			slave_mesh->all_second_order(false);

			// {
			// 	libMesh::MeshRefinement mesh_refinement(*master_mesh);
			// 	mesh_refinement.make_flags_parallel_consistent();
			// 	mesh_refinement.uniformly_refine(1);
			// }


			{
				libMesh::MeshRefinement mesh_refinement(*slave_mesh);
				mesh_refinement.make_flags_parallel_consistent();
				mesh_refinement.uniformly_refine(1);
			}

		} else {

			if(two_d) {
				auto elem_type_master  = libMesh::TRI6;
				auto elem_type_slave   = libMesh::QUAD8;

				libMesh::MeshTools::Generation::build_square(
					*master_mesh,
					n_master, n_master,
					0, 1.,
					0, 1.,
					elem_type_master
					);

				libMesh::MeshTools::Generation::build_square(
					*slave_mesh,
					n_slave, n_slave,
					0, 1.,
					0, 1.,
					elem_type_slave
					);
			} else {
				auto elem_type_master  = libMesh::TET10;
				auto elem_type_slave   = libMesh::TET10;

				libMesh::MeshTools::Generation::build_cube(
					*master_mesh,
					n_master, n_master, n_master,
					0, 1.,
					0, 1.,
					0, 1.,
					elem_type_master
					);

				libMesh::MeshTools::Generation::build_cube(
					*slave_mesh,
					n_slave, n_slave, n_slave,
					0, 1.,
					0, 1.,
					0, 1.,
					elem_type_slave
					);
			}
		}

		//equations system
		auto master_equation_systems = std::make_shared<libMesh::EquationSystems>(*master_mesh);
		auto &master_sys = master_equation_systems->add_system<libMesh::LinearImplicitSystem>("master_sys");

		auto slave_equation_systems = std::make_shared<libMesh::EquationSystems>(*slave_mesh);
		auto &slave_sys = slave_equation_systems->add_system<libMesh::LinearImplicitSystem>("slave_sys");

		//scalar function space
		auto V_m = FunctionSpaceT(master_equation_systems, libMesh::LAGRANGE, elem_order_master, "u_m");
		auto V_s = FunctionSpaceT(slave_equation_systems, libMesh::LAGRANGE,  elem_order_slave,  "u_s");

		V_m.initialize();
		V_s.initialize();

		Chrono c;
		c.start();
		USparseMatrix B;
		moonolith::Communicator comm(init.comm().get());

		const bool use_interpolation = true;
		if(assemble_volume_transfer(
			comm,
			master_mesh,
			slave_mesh,
			make_ref(V_m.dof_map()),
			make_ref(V_s.dof_map()),
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

			USparseMatrix T;
			if(use_interpolation) {
				T = B;
				UVector t = sum(T, 1);
				double t_max = max(t);
				double t_min = min(t);

				assert(t_min >= -1e-8);
				assert(t_max <= (1 + 1e-8));
				// std::cout << "[" << t_min << ", " << t_max << "] subset of [0, 1]" << std::endl;
				Interpolator(make_ref(T)).describe(std::cout);

			} else {
				USparseMatrix D_inv = diag(1./sum(B, 1));
				T = D_inv * B;
			}

			UVector v_m = local_values(V_m.dof_map().n_local_dofs(), 1.);

			auto f_rhs = ctx_fun< std::vector<double> >([](const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
				using std::sin;

				const auto &pts = ctx.fe()[0]->get_xyz();

				const auto n = pts.size();
				std::vector<double> ret(n);

				for(std::size_t i = 0; i != n; ++i) {
					double x = pts[i](0) - 0.5;
					double y = pts[i](1) - 0.5;
					double z = pts[i](2);
					ret[i] = sin(x) * (x*x + y*y + z*z);
				}

				return ret;
			});

			auto u = trial(V_m);
			auto v = test(V_m);
			auto p_form = inner(f_rhs, v) * dX;
			auto m_form = inner(u, v) * dX;

			UVector scaled_sol;
			USparseMatrix mass_mat;

			utopia::assemble(p_form, scaled_sol);
			utopia::assemble(m_form, mass_mat);

			if(elem_order_master == libMesh::FIRST) {
				v_m = e_mul(1./sum(mass_mat, 1), scaled_sol);
			} else {
				Factorization<USparseMatrix, UVector>().solve(mass_mat, scaled_sol, v_m);
			}

			// v_m.set(1.);

			UVector v_s = T * v_m;

			convert(v_s, *slave_sys.solution);
			slave_sys.solution->close();

			libMesh::Nemesis_IO surf_IO(*slave_mesh);
			surf_IO.write_equation_systems("interpolation_slave.e", *slave_equation_systems);

			convert(v_m, *master_sys.solution);
			master_sys.solution->close();

			libMesh::Nemesis_IO vol_IO(*master_mesh);
			vol_IO.write_equation_systems("interpolation_master.e", *master_equation_systems);

		} else {
			assert(false);
		}
	}
}
