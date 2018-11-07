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
#include "utopia_Grid2MeshSurfaceTransferAssembler.hpp"

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

		        is_shell = mesh_.mesh().spatial_dimension() > mesh_.mesh().mesh_dimension();


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
		bool is_shell;
	};

	class OwnershipRanges {
	public:
		OwnershipRanges(const int rank, const int comm_size)
		: rank_(rank)
		{
			ownership_ranges_.resize(comm_size + 1, 0);
		}

		void init(const std::size_t n)
		{
			//domain decomposition
			std::size_t comm_size = ownership_ranges_.size() - 1;
			std::size_t n_local = n / comm_size;
			std::size_t n_remainder = n % comm_size;

			ownership_ranges_[1] = n_local + n_remainder;
			for(std::size_t i = 2; i <= comm_size; ++i) {
				ownership_ranges_[i] = n_local + ownership_ranges_[i - 1];
			}
		}

		inline long local_nodes_begin() const
		{
			return ownership_ranges_[rank_];
		}

		inline long local_nodes_end() const
		{
			return ownership_ranges_[rank_ + 1];
		}

		inline long n_local_dofs() const
		{
			return local_nodes_end() - local_nodes_begin();
		}

		const std::vector<long> &get()
		{
			return ownership_ranges_;
		}

	private:
		int rank_;
		std::vector<long> ownership_ranges_;
	};

	class Grid2 : public Grid<2>, public Configurable {
	public:
		Grid2(const int rank, const int comm_size)
		: ownership_ranges(rank, comm_size)
		{}

		void read(Input &is) override
		{
		    try {

		    	Vector box_min(0., 0.);
		    	Vector box_max(1., 1.);

		    	is.get("grid", [&](Input &is) {
			    	is.get("n-x", this->dims[0]);
			    	is.get("n-y", this->dims[1]);

			    	is.get("min-x", box_min[0]);
			    	is.get("min-y", box_min[1]);

			    	is.get("max-x", box_max[0]);
			    	is.get("max-y", box_max[1]);
		    	});

		    	Vector box_range = box_max - box_min;

		    	this->map = [=](const Vector &x) -> Vector {
		    		//parametrization from unit-square to transformed quadrilateral
		    		Vector y = x;
		    		y.x *= box_range.x;
		    		y.y *= box_range.y;
		    		return y + box_min; 
		    	};

		    	this->inverse_map = [=](const Vector &x) -> Vector {
		    		//inverse-parametrization from  transformed quadrilateral to unit-square
		    		Vector y = x - box_min;
		    		y.x /= box_range.x;
		    		y.y /= box_range.y;
		    		return y;
		    	};

		    	this->dof_map = [](const Integer &hash) -> Integer {
		    		//dof-map starting from the hash (generated from x major)
		    		return hash;
		    	};

		    	ownership_ranges.init(this->n_nodes());

		    } catch(const std::exception &ex) {
		        std::cerr << ex.what() << std::endl;
		        assert(false);
		    }
		}

		inline long local_nodes_begin() const
		{
			return ownership_ranges.local_nodes_begin();
		}

		inline long local_nodes_end() const
		{
			return ownership_ranges.local_nodes_end();
		}

		inline long n_local_dofs() const
		{
			return ownership_ranges.n_local_dofs();
		}

		inline bool empty() const
		{
			return this->n_nodes() == 0;
		}

		OwnershipRanges ownership_ranges;
	};

	class Grid3 : public Grid<3>, public Configurable {
	public:
		Grid3(const int rank, const int comm_size)
		: ownership_ranges(rank, comm_size)
		{}

		void read(Input &is) override
		{
		    try {

		    	Vector box_min(0., 0., 0.);
		    	Vector box_max(1., 1., 1.);

		    	is.get("grid", [&](Input &is) {
			    	is.get("n-x", this->dims[0]);
			    	is.get("n-y", this->dims[1]);
			    	is.get("n-z", this->dims[2]);

			    	is.get("min-x", box_min[0]);
			    	is.get("min-y", box_min[1]);
			    	is.get("min-z", box_min[2]);

			    	is.get("max-x", box_max[0]);
			    	is.get("max-y", box_max[1]);
			    	is.get("max-z", box_max[2]);
		    	});

		    	Vector box_range = box_max - box_min;

		    	this->map = [=](const Vector &x) -> Vector {
		    		//parametrization from unit-cube to transformed hexahedron
		    		Vector y = x;
		    		y.x *= box_range.x;
		    		y.y *= box_range.y;
		    		y.z *= box_range.z;
		    		return y + box_min; 
		    	};

		    	this->inverse_map = [=](const Vector &x) -> Vector {
		    		//inverse-parametrization from  transformed hexahedron to unit-square
		    		Vector y = x - box_min;
		    		y.x /= box_range.x;
		    		y.y /= box_range.y;
		    		y.z /= box_range.z;
		    		return y;
		    	};

		    	this->dof_map = [](const Integer &hash) -> Integer {
		    		//dof-map starting from the hash (generated from x major)
		    		return hash;
		    	};

		    	ownership_ranges.init(this->n_nodes());

		    } catch(const std::exception &ex) {
		        std::cerr << ex.what() << std::endl;
		        assert(false);
		    }
		}

		inline long local_nodes_begin() const
		{
			return ownership_ranges.local_nodes_begin();
		}

		inline long local_nodes_end() const
		{
			return ownership_ranges.local_nodes_end();
		}

		inline long n_local_dofs() const
		{
			return ownership_ranges.n_local_dofs();
		}

		inline bool empty() const
		{
			return this->n_nodes() == 0;
		}

		OwnershipRanges ownership_ranges;
	};

	class Grid2MeshTransferApp::InputGrid : public Configurable {
	public:
		InputGrid(libMesh::Parallel::Communicator &comm)
		: comm_(comm),
		  grid2(comm.rank(), comm.size()),
		  grid3(comm.rank(), comm.size())
		{}

		void read(Input &is) override
		{
		    try {
		    	int dim = 3;
		    	is.get("dim", dim);

		    	if(dim == 3) {
		    		is_3D = true;
		    		grid3.read(is);
		    	} else {
		    		is_3D = false;
		    		assert(dim == 2);
		    		grid2.read(is);
		    	}

		    } catch(const std::exception &ex) {
		        std::cerr << ex.what() << std::endl;
		        assert(false);
		    }
		}

		libMesh::Parallel::Communicator &comm_;

		bool is_3D;
		Grid2 grid2;
		Grid3 grid3;
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
		bool enumerate_nodes = false;
		std::shared_ptr<LocalAssembler> local_assembler_;
		std::shared_ptr<Local2Global> local2global_;
		std::shared_ptr<TransferOperator> transfer_op_;

#ifdef WITH_TINY_EXPR
		std::shared_ptr<SymbolicFunction> fun;
#else
		std::shared_ptr<ConstantCoefficient<double, 0>> fun;
#endif //WITH_TINY_EXPR

		bool fun_is_constant;
		bool print_info = false;
		bool volume_to_surface = false;

		auto master_mesh_dim = 3;
		auto master_spatial_dim = 3;

		is_ptr->get("transfer", [&](Input &is) {
			//get spaces
			is_ptr->get("master", input_master);
			is_ptr->get("slave",  input_slave);

			if(!input_master.is_3D) {
				master_mesh_dim = 2;
				master_spatial_dim = 2;
			}

			is.get("write-operators-to-disk", write_operators_to_disk);
			is.get("type", type);
			is.get("assemble-mass-mat", assemble_mass_mat_);
			is.get("print-info", print_info);
			is.get("volume-to-surface", volume_to_surface);

			if(type == "l2-projection") {
				biorth_basis = true;
				is.get("biorth-basis", biorth_basis);
				
				local_assembler_ = std::make_shared<L2LocalAssembler>(
					master_mesh_dim,
					biorth_basis,
					assemble_mass_mat_ || !biorth_basis,
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
				local_assembler_ = apl2;
			}

			if(!local_assembler_) {
				assert(false);
				std::cerr << "choose type of assembler" << std::endl;
				return;
			}

			local2global_ = std::make_shared<Local2Global>(is_interpolation_);

			is.get("enumerate-nodes", enumerate_nodes);


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
		
		
		bool ok = false;

		if(input_master.is_3D) {
			auto &grid3 = input_master.grid3;

			if(print_info) {
				std::cout << "3D grid" << std::endl;
				grid3.describe(std::cout);
			}

			if(volume_to_surface) {
				Grid2MeshSurfaceTransferAssembler transfer_assembler(local_assembler_, local2global_);
				ok = transfer_assembler.assemble(
					grid3,
					grid3.ownership_ranges.get(),
					make_ref(input_slave.mesh()),
					make_ref(input_slave.space().dof_map()),
					mats,
					opts
				);

			} else {
				Grid2MeshTransferAssembler transfer_assembler(local_assembler_, local2global_);
				ok = transfer_assembler.assemble(
					grid3,
					grid3.ownership_ranges.get(),
					make_ref(input_slave.mesh()),
					make_ref(input_slave.space().dof_map()),
					mats,
					opts
				);
			}
		} else {
			auto &grid2 = input_master.grid2;

			if(print_info) {
				std::cout << "2D grid" << std::endl;
				grid2.describe(std::cout);
			}

			if(volume_to_surface) {
				Grid2MeshSurfaceTransferAssembler transfer_assembler(local_assembler_, local2global_);
				ok = transfer_assembler.assemble(
					grid2,
					grid2.ownership_ranges.get(),
					make_ref(input_slave.mesh()),
					make_ref(input_slave.space().dof_map()),
					mats,
					opts
				);
			} else {
				Grid2MeshTransferAssembler transfer_assembler(local_assembler_, local2global_);
				ok = transfer_assembler.assemble(
					grid2,
					grid2.ownership_ranges.get(),
					make_ref(input_slave.mesh()),
					make_ref(input_slave.space().dof_map()),
					mats,
					opts
				);
			}
		}

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
		

		OwnershipRanges * ownr = nullptr;

		if(input_master.is_3D) {
			ownr = &input_master.grid3.ownership_ranges;
		} else {
			ownr = &input_master.grid2.ownership_ranges;
		}

		fun_master = local_zeros(ownr->n_local_dofs());

		{
			Write<UVector> w_(fun_master);

			for(int i = ownr->local_nodes_begin(); i < ownr->local_nodes_end(); ++i) {
				if(enumerate_nodes) {
					fun_master.set(i, comm_->rank());
				} else {
#ifdef WITH_TINY_EXPR
					if(input_master.is_3D) {
						auto p = input_master.grid3.point(i);
						fun_master.set(i, fun->eval(p.x, p.y, p.z));
					} else {
						auto p = input_master.grid2.point(i);
						fun_master.set(i, fun->eval(p.x, p.y));
					}
#else
					fun_master.set(i, fun->expr());
#endif //WITH_TINY_EXPR
				}
			}
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

