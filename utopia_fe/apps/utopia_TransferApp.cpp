#include "utopia_TransferApp.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "libmesh/mesh_refinement.h"

namespace utopia {

	static void refine(const int n_refs, libMesh::MeshBase &mesh)
	{
		if(n_refs <= 0) return;

		libMesh::MeshRefinement mesh_refinement(mesh);
		mesh_refinement.make_flags_parallel_consistent();
		mesh_refinement.uniformly_refine(n_refs);
	}


	void TransferApp::init(libMesh::LibMeshInit &init)
	{
		comm_ = make_ref(init.comm());
	}

	void TransferApp::run(const std::string &conf_file_path)
	{
		Chrono c;

		c.start();

		mesh_master_ = std::make_shared<libMesh::DistributedMesh>(*comm_);
		mesh_slave_  = std::make_shared<libMesh::DistributedMesh>(*comm_);

		// std::cout << "running: " << conf_file_path << std::endl;
		auto is_ptr = open_istream(conf_file_path);

		is_ptr->read("transfer", [this](Input &is) {
			std::string path;
			type = "l2-projection"; //interpolation, approx-l2-projection
			int order = 1;
			std::string fe_family = "LAGRANGE";

			////////////////// MASTER ///////////////////////

			is.read("mesh-master", path);
			is.read("order-master", order);
			is.read("fe-family-master", fe_family);


			mesh_master_->read(path);

			if(order == 2) {
				mesh_master_->all_second_order(false);
			}

			int n_master_ref = 0;
			is.read("refine-master", n_master_ref);
			refine(n_master_ref, *mesh_master_);

			equation_systems_master_ = std::make_shared<libMesh::EquationSystems>(*mesh_master_);
			equation_systems_master_->add_system<libMesh::LinearImplicitSystem>("master");
			space_master_            = std::make_shared<LibMeshFunctionSpace>(equation_systems_master_, libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family), libMesh::Order(order), "u_master");
			space_master_->initialize();

			////////////////// SLAVE ///////////////////////
			order = 1;
			fe_family = "LAGRANGE";

			is.read("mesh-slave", path);
			is.read("order-slave", order);
			is.read("fe-family-slave", fe_family);
			is.read("type", type);
			mesh_slave_->read(path);

			if(order == 2) {
				mesh_slave_->all_second_order(false);
			}

			int n_slave_ref = 0;
			is.read("refine-slave", n_slave_ref);
			refine(n_slave_ref, *mesh_slave_);

			equation_systems_slave_ = std::make_shared<libMesh::EquationSystems>(*mesh_slave_);
			equation_systems_slave_->add_system<libMesh::LinearImplicitSystem>("slave");
			space_slave_            = std::make_shared<LibMeshFunctionSpace>(equation_systems_slave_, libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family), libMesh::Order(order), "u_slave");
			space_slave_->initialize();

			is_interpolation_ = false;
			assemble_mass_mat_ = 0;

			bool force_shell = false;
			is.read("force-shell", force_shell);

			is.read("assemble-mass-mat", assemble_mass_mat_);

			if(type == "l2-projection") {
				biorth_basis = true;
				is.read("biorth-basis", biorth_basis);
				
				local_assembler_ = std::make_shared<L2LocalAssembler>(
					mesh_master_->mesh_dimension(),
					biorth_basis,
					assemble_mass_mat_,
					force_shell || mesh_master_->mesh_dimension() < mesh_master_->spatial_dimension()
				);

			} else if(type == "interpolation") {
				local_assembler_ = std::make_shared<InterpolationLocalAssembler>(mesh_master_->mesh_dimension());
				is_interpolation_ = true;
			} else if(type == "approx-l2-projection") {
				int quad_order = -1;

				is.read("quad-order-approx", quad_order);
				std::cout << "quad_order: " << quad_order << std::endl;

				auto apl2 = std::make_shared<ApproxL2LocalAssembler>(mesh_master_->mesh_dimension());
				apl2->set_quadrature_order(quad_order);
				local_assembler_  = apl2;
			}

			if(!local_assembler_) {
				assert(false);
				std::cerr << "choose type of assembler" << std::endl;
				return;
			}

			local2global_ = std::make_shared<Local2Global>(is_interpolation_);

#ifdef WITH_TINY_EXPR
			std::string expr = "x";
			is.read("function", expr);

			fun = std::make_shared<SymbolicFunction>(expr);

			std::string fun_type = "non-linear";
			is.read("function-type", fun_type);

			if(fun_type == "constant") {
				fun_is_constant = true;
			} else {
				fun_is_constant = false;
			}
#else
			double expr = 1.;
			is.read("function", expr);
			fun_is_constant = true;

			fun = std::make_shared<ConstantCoefficient<double, 0>>(expr);
#endif //WITH_TINY_EXPR

		});

		c.stop();

		if(mpi_world_rank() == 0) {
			std::cout << "set-up time: " << c << std::endl;
		}


		TransferOptions opts;
		opts.from_var_num = 0;
		opts.to_var_num   = 0;
		opts.n_var        = 1;
		opts.tags         = {};

		//////////////////////////////////////////////////////

		c.start();

		std::vector<std::shared_ptr<USparseMatrix>> mats;
		TransferAssembler transfer_assembler(local_assembler_, local2global_);
		bool ok = transfer_assembler.assemble(
			mesh_master_,
			make_ref(space_master_->dof_map()),
			mesh_slave_,
			make_ref(space_slave_->dof_map()),
			mats,
			opts
		);

		if(!ok) {
			std::cerr << "[Error] transfer failed" << std::endl;
			return;
		}

		c.stop();

		if(mpi_world_rank() == 0) {
			std::cout << "assembly time: " << c << std::endl;
		}

		for(auto mat_ptr : mats) {
			double sum_m = sum(*mat_ptr);
			if(mpi_world_rank() == 0) {
				std::cout << "rows x cols = " << size(*mat_ptr).get(0) << " x " << size(*mat_ptr).get(1) << std::endl;
				std::cout << "sum(M): " << sum_m << std::endl;
			}
		}

		if(type == "l2-projection" || type == "approx-l2-projection") {
			if(biorth_basis && type != "approx-l2-projection") {
				auto pl2 = std::make_shared<PseudoL2TransferOperator>();
				pl2->init_from_coupling_operator(*mats[0]);
				transfer_op_ = pl2;
			} else {
				if(mats.size() == 2) {
					auto l2op = std::make_shared<L2TransferOperator>(mats[0], mats[1], std::make_shared<Factorization<USparseMatrix, UVector>>());
					l2op->fix_mass_matrix_operator();
					transfer_op_ = l2op;
				} else {
					auto u = trial(*space_slave_);
					auto v = test(*space_slave_);

					auto D = std::make_shared<USparseMatrix>();

					assemble(inner(u, v) * dX, *D);
					transfer_op_ = std::make_shared<L2TransferOperator>(mats[0], D, std::make_shared<Factorization<USparseMatrix, UVector>>());
				}
			}

		} else if(type == "interpolation") {
			transfer_op_ = std::make_shared<Interpolator>(mats[0]);
		}

		////////////////////////////////////////////////////////////

		auto u = trial(*space_master_);
		auto v = test(*space_master_);

		c.start();

		UVector fun_master_h, fun_master, fun_slave, back_fun_master;

		USparseMatrix mass_mat_master;
		assemble(inner(u, v) * dX, mass_mat_master);

		if(!fun_is_constant) {
			assemble(inner(*fun, v) * dX, fun_master_h);

			fun_master = fun_master_h;

			BiCGStab<USparseMatrix, UVector> solver;
			solver.solve(mass_mat_master, fun_master_h, fun_master);
		} else {
#ifdef WITH_TINY_EXPR
			fun_master = local_values(space_master_->dof_map().n_local_dofs(), fun->eval(0., 0., 0.));
#else
			fun_master = local_values(space_master_->dof_map().n_local_dofs(), fun->expr());
#endif //WITH_TINY_EXPR
		}

		c.stop();

		if(mpi_world_rank() == 0) {
			std::cout << "Assembled M and fun_m" << std::endl;
			std::cout << c << std::endl;
		}

		transfer_op_->apply(fun_master, fun_slave);
		transfer_op_->apply_transpose(fun_slave, back_fun_master);

		double sum_fun_master = sum(fun_master);
		double sum_fun_slave  = sum(fun_slave);

		transfer_op_->describe(std::cout);
		std::cout << "f_m = " << sum_fun_master << ", f_s = " << sum_fun_slave << std::endl;
		std::cout << "sum(M_m) = " << double(sum(mass_mat_master)) << std::endl;

		////////////////////////////////////////////////////////////
		//output
		////////////////////////////////////////////////////////////

		convert(fun_master, *space_master_->equation_system().solution);
		space_master_->equation_system().solution->close();

		libMesh::Nemesis_IO io_master(*mesh_master_);
		io_master.write_equation_systems("master.e", *equation_systems_master_);

		convert(back_fun_master, *space_master_->equation_system().solution);
		space_master_->equation_system().solution->close();

		libMesh::Nemesis_IO io_master_adj(*mesh_master_);
		io_master_adj.write_equation_systems("master_adj.e", *equation_systems_master_);

		////////////////////////////////////////////////////////////
		convert(fun_slave, *space_slave_->equation_system().solution);
		space_slave_->equation_system().solution->close();
		libMesh::Nemesis_IO io_slave(*mesh_slave_);
		io_slave.write_equation_systems("slave.e", *equation_systems_slave_);
	}

	TransferApp::TransferApp() {}
	TransferApp::~TransferApp() {}
}

