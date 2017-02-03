#include "utopia_FEDSLMortarExamples.hpp"

#include <iostream>
#include "utopia.hpp"

//fe extension
#include "utopia_fe.hpp"
#include "MortarAssembler.hpp"
#include "MixedParMortarAssembler.hpp"
#include "ParMortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>

#include "utopia_ContactSimParams.hpp"
#include "libmesh/linear_partitioner.h"


#include "express_Profiler.hpp"
#include "utopia_Polygon.hpp"

using namespace utopia;
using namespace std;
using namespace libMesh;


namespace utopia {


	void mortar_transfer_aux(const std::shared_ptr<Mesh> &mesh_master, 
							 const std::shared_ptr<Mesh> &mesh_slave,
							 const libMesh::Order order_elem = FIRST)
	{
		Chrono c;
		c.start();

		
		
		int order_quad = order_elem + order_elem + order_elem;

		LibMeshFEContext<LinearImplicitSystem> master_context(mesh_master);
		auto master_space = fe_space(LAGRANGE, order_elem, master_context);
		master_context.equation_systems.init();

		LibMeshFEContext<LinearImplicitSystem> slave_context(mesh_slave);
		auto slave_space = fe_space(LAGRANGE, order_elem, slave_context);
		slave_context.equation_systems.init();

		auto v = make_ref(master_space);
		bool use_biorthogonal_mults = true;

		//////////////////////////////////////////////////
		//////////////////////////////////////////////////


		MortarAssembler assembler(make_ref(master_space), make_ref(slave_space));
		assembler.set_use_biorthogonal_multipliers(use_biorthogonal_mults);


		DSMatrixd matrix;
		assembler.assemble(matrix);

		//////////////////////////////////////////////////
		//////////////////////////////////////////////////
		
		int dim = mesh_slave->mesh_dimension();
		auto u = fe_function(slave_space);
		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, libMesh::Order(order_quad)));

		Size s = size(matrix);
		DSMatrixd mass = sparse(s.get(0), s.get(0), 8);
		{
			Write<DSMatrixd> w_mass(mass);
			assemble(u, u, integral(dot(u, u)), mass);
		}

		DVectord master_fun   = values(s.get(1), 0.0);
		DSMatrixd master_mass = sparse(s.get(1), s.get(1), 8);
		{
			Write<DVectord> w_fun(master_fun);

			auto w = fe_function(master_space);
			w.set_quad_rule(make_shared<libMesh::QGauss>(dim, libMesh::Order(order_quad)));

			std::function<Real(const Point &p) > f = [dim](const Point &p) -> Real {
				Real center = 0;

				Real ret = 0;
				for(int i = 0; i < dim; ++i) {
					ret += (p(i) - center) * (p(i) - center);
				}

				return 10 * std::sqrt(ret) - 5;
			};

			assemble(w, integral(dot(coeff(f), w)), master_fun);

			Write<DSMatrixd> w_mm(master_mass);
			assemble(w, w, integral(dot(w, w)), master_mass);
		}




		// DVectord val_master = values(s.get(1), 1.0);
		DVectord val_master = zeros(s.get(1));
		solve(master_mass, master_fun, val_master);

		DVectord val_mult   = matrix * val_master;
		DVectord val_slave  = zeros(s.get(0), s.get(0));
		// disp(val_slave);
		// solve(mass, val_mult, val_slave);

		solve(assembler.D, val_mult, val_slave);

		{
			Read<DVectord> rl(val_slave);
			std::cout << "===\n" << val_slave.get(0) << std::endl;

		}

		// disp(assembler.D);


		c.stop();
		std::cout << "|master_elements| = " << mesh_master->n_elem() << std::endl;
		std::cout << "|slave_elements|  = " << mesh_slave->n_elem()  << std::endl;
		std::cout << double(sum(mass)) << "==" << double(sum(matrix)) << std::endl;
		std::cout << "time: ";
		c.describe(std::cout);

		// disp(val_slave);

		// DVectord m_sum = sum(matrix, 1);
		// disp(m_sum);

		// DSMatrixd T = diag(1./sum(matrix, 1)) * (matrix);
		// val_slave = T * val_master;

		std::cout << "no mass matrix\n--------------------------" << std::endl;
		// disp(val_slave);

		convert(val_slave, *slave_context.system.solution);
		ExodusII_IO(*mesh_slave).write_equation_systems ("transfer_2D_slave.e", slave_context.equation_systems);

		convert(val_master, *master_context.system.solution);
		ExodusII_IO(*mesh_master).write_equation_systems ("transfer_2D_master.e", master_context.equation_systems);

		// write("ser.m", matrix);
	}



    
    void mixed_par_mortar_transfer_aux(libMesh::Parallel::Communicator &libmesh_comm, const std::shared_ptr<Mesh> &mesh_master, const std::shared_ptr<Mesh> &mesh_slave)
    {
//         EXPRESS_EVENT_BEGIN("spaces");
//        
//         Chrono c;
//         c.start();
        
        auto order_elem = FIRST;
        int order_quad = order_elem + order_elem;
        
        LibMeshFEContext<LinearImplicitSystem> master_context(mesh_master);
        auto master_space = fe_space(LAGRANGE, order_elem, master_context);
        master_context.equation_systems.init();
        
        LibMeshFEContext<LinearImplicitSystem> slave_context(mesh_slave);
        auto slave_space = fe_space(LAGRANGE, order_elem, slave_context);
        slave_context.equation_systems.init();
        
        bool use_biorthogonal_mults = true;
        
       
        
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        
        //        make_ref(master_space)->system();
        //        std::cout<<"I am a master system"<<make_ref(master_space)->system()<<std::endl;
        
        // MPI_Comm comm = MPI_COMM_WORLD;
        MixedParMortarAssembler assembler(libmesh_comm, make_ref(master_space), make_ref(slave_space));
        assembler.set_use_biorthogonal_multipliers(use_biorthogonal_mults);
        
        // std::cout<<"I am a slave system"<<make_ref(slave_space)->system()<<std::endl;
        
        
        
        //		EXPRESS_EVENT_END("spaces");
        
        DSMatrixd transfer_op;
        DSMatrixd matrix;
        assembler.Assemble(matrix);
        assembler.Transfer(matrix,transfer_op);
        
        
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        
        //      EXPRESS_EVENT_BEGIN("mass_matrix");
        
        const Size gs = size(matrix);
        const Size ls = local_size(matrix);
        
        int dim = mesh_slave->mesh_dimension();
        auto u = fe_function(slave_space);
        u.set_quad_rule(make_shared<libMesh::QGauss>(dim, libMesh::Order(order_quad)));
        
        auto &mass_matrix = *slave_context.system.matrix;
        auto ass = make_assembly([&]() -> void {
         	assemble(u, u, integral(inner(u, u)), mass_matrix);
         });
        
        slave_context.system.attach_assemble_object(ass);
        slave_context.equation_systems.print_info();
        slave_context.equation_systems.solve();
     
        DSMatrixd mass;
        convert(mass_matrix, mass);
        
        DVectord val_master = local_values(ls.get(1), 1.0);
        DVectord val_mult   = matrix * val_master;
        DVectord val_slave  = local_zeros(ls.get(0));
        
        solve(mass, val_mult, val_slave);
        
        DVectord val_slave_pseudo=transfer_op*val_master;
        
//        monitor(0, matrix);
        
        // c.stop();
        std::cout << "|master_elements| = " << mesh_master->n_active_elem() << std::endl;
        std::cout << "|slave_elements|  = " << mesh_slave->n_active_elem()  << std::endl;
//        std::cout << double(sum(mass)) << "==" << double(sum(matrix)) << std::endl;
        // std::cout << "time: ";
        // c.describe(std::cout);
        
        
        //	EXPRESS_EVENT_END("mass_matrix");
        
        //disp(val_slave);
       // disp(val_slave_pseudo);
        // write("par.m", matrix);
    }
    
    
	void par_mortar_transfer_aux(libMesh::Parallel::Communicator &libmesh_comm,const std::shared_ptr<MeshBase> &master_slave)
	{
		Chrono c;
		c.start();

		auto order_elem = FIRST;
		int order_quad = order_elem + order_elem;

		LibMeshFEContext<LinearImplicitSystem> master_slave_context(master_slave);
		auto master_slave_space = fe_space(LAGRANGE, order_elem, master_slave_context);
		master_slave_context.equation_systems.init();

        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
		
		ParMortarAssembler assembler(libmesh_comm, make_ref(master_slave_space));

		DSMatrixd matrix;

		// EXPRESS_EVENT_BEGIN("l2assembly");
		assembler.Assemble(matrix);
		// EXPRESS_EVENT_END("l2assembly");
	}

    
    
    
    
    
    void par_mortar_surface_transfer_aux(libMesh::Parallel::Communicator &libmesh_comm,const std::shared_ptr<MeshBase> &master_slave)
    {
        Chrono c;
        c.start();
        
        auto order_elem = FIRST;
        int order_quad = order_elem + order_elem;
        
        LibMeshFEContext<LinearImplicitSystem> master_slave_context(master_slave);
        auto master_slave_space = fe_space(LAGRANGE, order_elem, master_slave_context);
        master_slave_context.equation_systems.init();
        
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        
        ParMortarAssembler surface_assembler(libmesh_comm, make_ref(master_slave_space));
        
        DSMatrixd matrix;
        const  libMesh::Real search_radius=0.4;
        // EXPRESS_EVENT_BEGIN("l2assembly");
        surface_assembler.SurfaceAssemble(matrix,search_radius,101,102);
        // EXPRESS_EVENT_END("l2assembly");
    }

    
    
	void mortar_transfer_2D_monolithic(LibMeshInit &init)
	{
		auto mesh = make_shared<Mesh>(init.comm());
        EXPRESS_EVENT_BEGIN("set_up");
		//mesh->partitioner().reset(new LinearPartitioner());
		mesh->read("../data/master_slave2D_new.e");
		par_mortar_transfer_aux(init.comm(),mesh);
        EXPRESS_EVENT_END("set_up");
	}

	void mortar_transfer_2D(LibMeshInit &init)
	{
		std::cout << "-----------------------------\n";
		std::cout << "mortar_transfer_2D\n";
		//////////////////////////////////////////////////
		//////////////////////////////////////////////////
		int n_master = 30; 
		int n_slave  = 10; 

		auto mesh_master = make_shared<Mesh>(init.comm());		
		MeshTools::Generation::build_square (*mesh_master,
			n_master, n_master,
			-1, 1.,
			-1, 1.,
			QUAD4);

		//////////////////////////////////////////////////
		//////////////////////////////////////////////////

		auto mesh_slave = make_shared<Mesh>(init.comm());		
		MeshTools::Generation::build_square (*mesh_slave,
			n_slave, n_slave,
			-1, 1.,
			-1, 1.,
			TRI3);

		// MeshTools::Modification::smooth (*mesh_slave,
		//                                        10,
		//                                        0.4);

		// MeshTools::Modification::distort (*mesh_slave,
		//                                   0.03,
		//                                   true);

		// MeshTools::Modification::scale (*mesh_master,
		//   1.1,
		//   1.1,
		//   1.1);
		  
		//   MeshTools::Modification::distort (*mesh_master,
		//                                     0.01,
		//                                     true);
		 
		// plot_polygon_mesh(*mesh_master, "mater");
		// plot_polygon_mesh(*mesh_slave, "slave");


//		plot_mesh(*mesh_master, "mater");
//		plot_mesh(*mesh_slave, "slave");


		mixed_par_mortar_transfer_aux(init.comm(), mesh_master, mesh_slave);
		//mortar_transfer_aux(mesh_master, mesh_slave);
		std::cout << "-----------------------------\n";

	}

	void mortar_transfer_3D_monolithic(LibMeshInit &init)
	{

		EXPRESS_EVENT_BEGIN("set_up");
		auto mesh = make_shared<Mesh>(init.comm());

		//mesh->partitioner().reset(new LinearPartitioner());
            // Read the mesh file. Here the file lshape.unv contains
            // an L--shaped domain in .unv format.
       // mesh->read("../data/cube12_space5.e"); //("../data/master_slave3D_translated.e");
       //mesh->read("../data/master_slave3D_translated.e");
        mesh->read("../data/rect.e");

            // Print information about the mesh to the screen.
		// mesh->print_info();

		EXPRESS_EVENT_END("set_up");

		//par_mortar_transfer_aux(init.comm(),mesh);
        par_mortar_surface_transfer_aux(init.comm(),mesh);

		std::cout << "-----------------------------\n";
	}

    void mortar_transfer_3D(LibMeshInit &init)
    {
        //EXPRESS_EVENT_BEGIN("set_up");
        std::cout << "-----------------------------\n";
        std::cout << "mortar_transfer_3D\n";
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        int n_master = 2;
        int n_slave  = 3;
        
        auto mesh_master = make_shared<Mesh>(init.comm());
        
        //mesh_master->partitioner().reset(new SFCPartitioner());
;
        
        
        MeshTools::Generation::build_cube(*mesh_master,
                                          n_master, n_master, n_master,
                                          -2., 3.,
                                          -2., 3.,
                                          -2., 3.,
                                          TET4);
        
        
        
        //////////////////////////////////////////////////
        //////////////////////////////////////////////////
        
        
        
        auto mesh_slave = make_shared<Mesh>(init.comm());
        
        //mesh_slave->partitioner().reset(new SFCPartitioner());
        
        MeshTools::Generation::build_cube (*mesh_slave,
                                           n_slave, n_slave, n_slave,
                                           -2., 3.,
                                           -2., 3.,
                                           -2., 3.,
                                           TET4);
        
        
        
        //EXPRESS_EVENT_END("set_up");
        //mortar_transfer_aux(mesh_master, mesh_slave);
        mixed_par_mortar_transfer_aux(init.comm(), mesh_master, mesh_slave);
        
        
        //std::cout << "-----------------------------\n";
    }

	void surface_mortar(LibMeshInit &init)
	{
		std::cout << "-----------------------------\n";
		std::cout << "surface_mortar\n";
		//////////////////////////////////////////////////
		//////////////////////////////////////////////////

		//ContactSimParams params = leaflets_contact;
		 ContactSimParams params = contact8;
    	// ContactSimParams params = multi_contact_quads;
		// ContactSimParams params = triple_contact_circle;
		// ContactSimParams params = contact_3D_tets;

		auto mesh = make_shared<Mesh>(init.comm());		
		mesh->read(params.mesh_path);
//		plot_mesh(*mesh, "mesh");
        


		const int dim = mesh->mesh_dimension();

		//////////////////////////////////////////////////
		//////////////////////////////////////////////////

		Chrono c;
		c.start();

		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		auto space = vector_fe_space("disp_", LAGRANGE_VEC, FIRST, context);
		auto u = fe_function(space);


		static const bool is_leaflet = true;

		//boundary conditions
		std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_1 = [dim, &params](const Point &, DenseVector<Real> &output) {
			output.resize(dim);
			output.zero();

			if(is_leaflet) {
				// output(0) = -params.dirichlet_value_1;
				output(1) = params.dirichlet_value_1;
			} else {
				//offset to y coordinate of displacement
				output(1) = params.dirichlet_value_1;
			}
		};

		auto bc_1 = boundary_conditions(u == vec_coeff(boundary_cond_1),    { params.boundary_tag_1 });
		strong_enforce(bc_1);

		std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_2 = [dim, &params](const Point &, DenseVector<Real> &output) {
			output.resize(dim);
			output.zero();

			if(is_leaflet) {
				output(0) = -params.dirichlet_value_2;
				output(1) = -params.dirichlet_value_2;
			} else {
				//offset to y coordinate of displacement
				output(1) = params.dirichlet_value_2;
			}
		};
		
		auto bc_2 = boundary_conditions(u == vec_coeff(boundary_cond_2), { params.boundary_tag_2 });
		strong_enforce(bc_2);



		std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_3 = [dim, &params](const Point &, DenseVector<Real> &output) {
			output.resize(dim);
			output.zero();

			if(is_leaflet) {
				output(0) = params.dirichlet_value_3;
				output(1) = -params.dirichlet_value_3;
			} else {
				//offset to y coordinate of displacement
				output(1) = params.dirichlet_value_3;
			}
		};

		if(params.boundary_tag_3 >= 0) {
			auto bc_3 = boundary_conditions(u == vec_coeff(boundary_cond_3), { params.boundary_tag_3 });
			strong_enforce(bc_3);
		}


	

		context.equation_systems.init();
		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));

		double mu   = 1.0, lambda = 1.0;
		auto e      = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
		auto b_form = integral((mu * 0.5) * dot(e, e) + lambda * dot(div(u), div(u)));

		DenseVector<Real> vec(dim);
		vec.zero();

		//linear forms
		auto f 	    = vec_coeff(vec);
		auto l_form = integral(dot(f, u));
		
		//FIXME 
		auto ass = make_assembly([&]() -> void {
			assemble(u, u, b_form, l_form, 
				*context.system.matrix,  
				*context.system.rhs);
		});

		context.system.attach_assemble_object(ass);
		context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
		context.equation_systems.solve();

		c.stop();
		std::cout << "Matrix and rhs assembly: ";
		c.describe(std::cout); 
		std::cout << std::endl;

		long n = u.dof_map().n_dofs();

		MortarContactAssembler assembler(make_ref(space));
		assembler.set_strict_gap_policy(params.strict_search_radius);

		DSMatrixd coupling;
		DVectord gap;
		DVectord normals;
		DSMatrixd orhtogonal_trafos;

		std::vector<bool> is_contact_node;
		if(!assembler.assemble(coupling, 
			gap, 
			normals, 
			orhtogonal_trafos, 
			is_contact_node, 
			params.search_radius))
		{
			//Just set some values that do not change the original system
			coupling = identity(n, n);
			orhtogonal_trafos = identity(n, n);
			gap = values(n, 10000);
			normals = zeros(n);
			is_contact_node.resize(n);
			std::fill(is_contact_node.begin(), is_contact_node.end(), false);
		}

		DSMatrixd K  = sparse(n, n, 20);
		DVectord rhs = zeros(n);

		convert(*context.system.rhs, rhs);
		convert(*context.system.matrix, K);

		//Change of basis
		DVectord  sol_c = zeros(size(rhs));
		DVectord  rhs_c = transpose(orhtogonal_trafos) * transpose(coupling) * rhs;
		DSMatrixd K_c   = transpose(orhtogonal_trafos) * DSMatrixd(transpose(coupling) * K * coupling) * orhtogonal_trafos;
		DVectord  gap_c = transpose(coupling) * gap;
		apply_boundary_conditions(u, K_c, rhs_c);

		SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<Factorization<DSMatrixd, DVectord> >());
        // SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<ConjugateGradient<DSMatrixd, DVectord> >());
		
		newton.verbose(true);
		newton.solve(sol_c, K_c, rhs_c, gap_c);

		//Change back to original basis
		DVectord sol = coupling * (orhtogonal_trafos * sol_c);

		convert(sol, *context.system.solution);
		ExodusII_IO(*mesh).write_equation_systems ("surface_mortar.e", context.equation_systems);

		{
			Read<DVectord> w_g(gap);
			Read<DVectord> w_n(normals);

			for(size_t i = 0; i < is_contact_node.size(); ++i) {
				if(is_contact_node[i]) {
					const double len = gap.get((i/dim)*dim);
					context.system.solution->set(i, len * normals.get(i));
				} else {
					context.system.solution->set(i, 0);
				}
			}
		}

		ExodusII_IO(*mesh).write_equation_systems ("surface_mortar_cn.e", context.equation_systems);
	}


	void run_curved_poly_disc()
	{
		int order = 3;
		auto f = [order](const double * x, double *fx) -> void {
			fx[0] = x[0] * x[0] * x[0];
			fx[1] = x[1] * x[1] * x[1];
		};

		const double from[2] = { -1., 0. };
		const double to[2]   = {  0., 1. };

		std::vector<double> params_points, polyline;
		discretize_curve<2>(f, from, to, order, params_points, polyline);

		std::cout << "-----------------------\n";
		for(int i = 0; i < polyline.size(); i += 2) {
			std::cout << polyline[i] << ", " << polyline[i+1] << std::endl;
		}
		std::cout << "-----------------------\n";
	//	plot_polygon(2, polyline.size()/2, &polyline[0], "polyline");
	}


	void run_mortar_examples(libMesh::LibMeshInit &init)
	{

		EXPRESS_PROFILING_BEGIN()


		//mortar_transfer_2D(init);
        //mortar_transfer_3D(init);
		mortar_transfer_3D_monolithic(init);
		//surface_mortar(init);

		//run_curved_poly_disc();

		EXPRESS_PROFILING_END();
	}

}

