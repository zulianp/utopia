#include "utopia_Grid2MeshTransferApp.hpp"
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

#include "utopia_GridMeshTransfer.hpp"

namespace utopia {

	class Grid2MeshTransferApp::InputSpace : public Configurable {
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

	class Grid2MeshTransferApp::InputGrid : public Configurable {
	public:
		InputGrid(libMesh::Parallel::Communicator &comm)
		: comm_(comm)
		{
			ownership_ranges.resize(comm_.size() + 1, 0);
		}

		void read(Input &is) override
		{
		    try {
		    	is.get("nx", grid.dims[0]);
		    	is.get("ny", grid.dims[1]);
		    	is.get("nz", grid.dims[2]);

		    	//FIXME
		    	ownership_ranges[1] = grid.n_nodes();

		    } catch(const std::exception &ex) {
		        std::cerr << ex.what() << std::endl;
		        assert(false);
		    }
		}

		inline long local_nodes_begin() const
		{
			return ownership_ranges[comm_.rank()];
		}


		inline long local_nodes_end() const
		{
			return ownership_ranges[comm_.rank() + 1];
		}

		inline long n_local_dofs() const
		{
			return local_nodes_end() - local_nodes_begin();
		}

		libMesh::Parallel::Communicator &comm_;
		Grid<3> grid;
		std::vector<long> ownership_ranges;
	};

	void Grid2MeshTransferApp::init(libMesh::LibMeshInit &init)
	{
		comm_ = make_ref(init.comm());
	}

	void Grid2MeshTransferApp::run(const std::string &conf_file_path)
	{
		Chrono c;
		c.start();

		auto is_ptr = open_istream(conf_file_path);
		InputGrid input_master(*comm_);
		InputSpace input_slave(*comm_);

		//get operator props
		std::string path;
		std::string type = "l2-projection"; //interpolation, approx-l2-projection
		int order = 1;
		bool write_operators_to_disk = false;
		bool is_interpolation_ = false;
		bool assemble_mass_mat_ = 0;
		bool force_shell  = false;
		bool biorth_basis = false;
		std::shared_ptr<LocalAssembler> local_assembler_;
		std::shared_ptr<Local2Global> local2global_;
		std::shared_ptr<TransferOperator> transfer_op_;

#ifdef WITH_TINY_EXPR
		std::shared_ptr<SymbolicFunction> fun;
#else
		std::shared_ptr<ConstantCoefficient<double, 0>> fun;
#endif //WITH_TINY_EXPR

		bool fun_is_constant;


		auto master_mesh_dim = 3;
		auto master_spatial_dim = 3;

		is_ptr->get("transfer", [&](Input &is) {
			//get spaces
			is_ptr->get("master", input_master);
			is_ptr->get("slave",  input_slave);

		

			is.get("write-operators-to-disk", write_operators_to_disk);
			is.get("type", type);
			is.get("assemble-mass-mat", assemble_mass_mat_);

			if(type == "l2-projection") {
				biorth_basis = true;
				is.get("biorth-basis", biorth_basis);
				
				local_assembler_ = std::make_shared<L2LocalAssembler>(
					master_mesh_dim,
					biorth_basis,
					assemble_mass_mat_,
					force_shell || master_mesh_dim < master_spatial_dim
				);

			} else if(type == "interpolation") {
				local_assembler_ = std::make_shared<InterpolationLocalAssembler>(master_mesh_dim);
				is_interpolation_ = true;
			} else if(type == "approx-l2-projection") {
				int quad_order = -1;

				is.get("quad-order-approx", quad_order);
				std::cout << "quad_order: " << quad_order << std::endl;

				auto apl2 = std::make_shared<ApproxL2LocalAssembler>(master_mesh_dim);
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
		GridMeshTransfer<3> transfer_assembler(local_assembler_, local2global_);
		bool ok = transfer_assembler.assemble(
			input_master.grid,
			input_master.ownership_ranges,
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
					l2op->fix_mass_matrix_operator();
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
		c.start();

		UVector fun_master, fun_slave;
		fun_master = local_zeros(input_master.n_local_dofs());

		for(int i = input_master.local_nodes_begin(); i < input_master.local_nodes_end(); ++i) {
#ifdef WITH_TINY_EXPR
			auto p = input_master.grid.point(i);
			fun_master.set(i, fun->eval(p.x, p.y, p.z));
#else
			fun_master.set(i, fun->expr());
#endif //WITH_TINY_EXPR
		}

		c.stop();

		if(mpi_world_rank() == 0) {
			std::cout << "Assembled M and fun_m" << std::endl;
			std::cout << c << std::endl;
		}

		transfer_op_->apply(fun_master, fun_slave);

		double sum_fun_master = sum(fun_master);
		double sum_fun_slave  = sum(fun_slave);

		transfer_op_->describe(std::cout);
		std::cout << "f_m = " << sum_fun_master << ", f_s = " << sum_fun_slave << std::endl;

		////////////////////////////////////////////////////////////
		//output
		////////////////////////////////////////////////////////////
		convert(fun_slave, *input_slave.space().equation_system().solution);
		input_slave.space().equation_system().solution->close();
		libMesh::Nemesis_IO io_slave(input_slave.mesh());
		io_slave.write_equation_systems("slave.e", input_slave.space().equation_systems());
	}

	Grid2MeshTransferApp::Grid2MeshTransferApp() {}
	Grid2MeshTransferApp::~Grid2MeshTransferApp() {}
}

