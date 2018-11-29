#include "utopia_TransferApp.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_ui.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIScalarSampler.hpp"

#include "libmesh/mesh_refinement.h"

namespace utopia {

	class TransferApp::InputSpace : public Configurable {
	public:
		InputSpace(libMesh::Parallel::Communicator &comm)
		: mesh_(comm), space_(make_ref(mesh_))
		{}

		void read(Input &is) override
		{
		    try {
		        is.get("mesh", mesh_);
		        is.get("space", space_);


		    } catch(const std::exception &ex) {
		        std::cerr << ex.what() << std::endl;
		        assert(false);
		    }
		}

		inline bool empty() const
		{
		    return mesh_.empty();
		}

		inline libMesh::MeshBase &mesh()
		{
			return mesh_.mesh();
		}


		inline LibMeshFunctionSpace &space()
		{
			return space_.space()[0];
		}

		UIMesh<libMesh::DistributedMesh> mesh_;
		UIFunctionSpace<LibMeshFunctionSpace>  space_;
	};

	void TransferApp::init(libMesh::LibMeshInit &init)
	{
		comm_ = make_ref(init.comm());
	}

	void TransferApp::run(const std::string &conf_file_path)
	{
		Chrono c;
		c.start();

		auto is_ptr = open_istream(conf_file_path);
		InputSpace input_master(*comm_);
		InputSpace input_slave(*comm_);
		double tol = 1e-16;

		is_ptr->get("transfer", [&](Input &is) {
			//get spaces
			is.get("master", input_master);
			is.get("slave",  input_slave);

			//get operator props
			std::string path;
			type = "l2-projection"; //interpolation, approx-l2-projection
			int order = 1;
			std::string fe_family = "LAGRANGE";
			write_operators_to_disk = false;
			is_interpolation_ = false;
			assemble_mass_mat_ = 0;
			bool force_shell = false;

			is.get("write-operators-to-disk", write_operators_to_disk);
			is.get("type", type);
			is.get("force-shell", force_shell);
			is.get("assemble-mass-mat", assemble_mass_mat_);
			is.get("tol", tol);

			if(type == "l2-projection") {
				biorth_basis = true;
				is.get("biorth-basis", biorth_basis);
				
				local_assembler_ = std::make_shared<L2LocalAssembler>(
					input_master.mesh().mesh_dimension(),
					biorth_basis,
					assemble_mass_mat_,
					force_shell || input_master.mesh().mesh_dimension() < input_master.mesh().spatial_dimension()
				);

			} else if(type == "interpolation") {
				local_assembler_ = std::make_shared<InterpolationLocalAssembler>(input_master.mesh().mesh_dimension());
				is_interpolation_ = true;
			} else if(type == "approx-l2-projection") {
				int quad_order = -1;

				is.get("quad-order-approx", quad_order);
				std::cout << "quad_order: " << quad_order << std::endl;

				auto apl2 = std::make_shared<ApproxL2LocalAssembler>(input_master.mesh().mesh_dimension());
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
			is.get("function", expr);

			fun = std::make_shared<SymbolicFunction>(expr);

			std::string fun_type = "non-linear";
			is.get("function-type", fun_type);

			if(fun_type == "constant") {
				fun_is_constant = true;
			} else {
				fun_is_constant = false;
			}
#else
			double expr = 1.;
			is.get("function", expr);
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
			make_ref(input_master.mesh()),
			make_ref(input_master.space().dof_map()),
			make_ref(input_slave.mesh()),
			make_ref(input_slave.space().dof_map()),
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

		int op_num = 0;
		for(auto mat_ptr : mats) {
			double sum_m = sum(*mat_ptr);
			if(mpi_world_rank() == 0) {
				std::cout << "rows x cols = " << size(*mat_ptr).get(0) << " x " << size(*mat_ptr).get(1) << std::endl;
				std::cout << "sum(M): " << sum_m << std::endl;
			}

			if(write_operators_to_disk) {
				write("M" + std::to_string(op_num) + ".m", *mat_ptr);
			}

			++op_num;
		}

		if(type == "l2-projection" || type == "approx-l2-projection") {
			if(biorth_basis && type != "approx-l2-projection") {
				auto pl2 = std::make_shared<PseudoL2TransferOperator>();
				pl2->init_from_coupling_operator(*mats[0]);
				transfer_op_ = pl2;
			} else {
				if(mats.size() == 2) {
					auto l2op = std::make_shared<L2TransferOperator>(mats[0], mats[1], std::make_shared<Factorization<USparseMatrix, UVector>>());
					l2op->fix_mass_matrix_operator(tol);
					transfer_op_ = l2op;
				} else {
					auto u = trial(input_slave.space());
					auto v = test(input_slave.space());

					auto D = std::make_shared<USparseMatrix>();

					assemble(inner(u, v) * dX, *D);
					transfer_op_ = std::make_shared<L2TransferOperator>(mats[0], D, std::make_shared<Factorization<USparseMatrix, UVector>>());
				}
			}

		} else if(type == "interpolation") {
			transfer_op_ = std::make_shared<Interpolator>(mats[0]);
		}

		////////////////////////////////////////////////////////////

		auto u = trial(input_master.space());
		auto v = test(input_master.space());

		c.start();

		UVector fun_master_h, fun_master, fun_slave, back_fun_master;

		USparseMatrix mass_mat_master;
		assemble(inner(u, v) * dX, mass_mat_master);

		if(write_operators_to_disk) {
			write("M_m.m", mass_mat_master);
		}

		if(!fun_is_constant) {
			assemble(inner(*fun, v) * dX, fun_master_h);

			fun_master = fun_master_h;

			BiCGStab<USparseMatrix, UVector> solver;
			solver.solve(mass_mat_master, fun_master_h, fun_master);
		} else {
#ifdef WITH_TINY_EXPR
			fun_master = local_values(input_master.space().dof_map().n_local_dofs(), fun->eval(0., 0., 0.));
#else
			fun_master = local_values(input_master.space().dof_map().n_local_dofs(), fun->expr());
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

		convert(fun_master, *input_master.space().equation_system().solution);
		input_master.space().equation_system().solution->close();

		libMesh::Nemesis_IO io_master(input_master.mesh());
		io_master.write_equation_systems("master.e", input_master.space().equation_systems());

		convert(back_fun_master, *input_master.space().equation_system().solution);
		input_master.space().equation_system().solution->close();

		libMesh::Nemesis_IO io_master_adj(input_master.mesh());
		io_master_adj.write_equation_systems("master_adj.e", input_master.space().equation_systems());

		////////////////////////////////////////////////////////////
		convert(fun_slave, *input_slave.space().equation_system().solution);
		input_slave.space().equation_system().solution->close();
		libMesh::Nemesis_IO io_slave(input_slave.mesh());
		io_slave.write_equation_systems("slave.e", input_slave.space().equation_systems());
	}

	TransferApp::TransferApp() {}
	TransferApp::~TransferApp() {}
}

