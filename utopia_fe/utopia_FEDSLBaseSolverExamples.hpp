#ifndef FEDSL_SOLVER_EXAMPLES_HPP
#define FEDSL_SOLVER_EXAMPLES_HPP



#include <iostream>
#include <cmath>
#include "utopia_FEDSLBaseSolverExamples.hpp"
#include "utopia.hpp"

//fe extension
#include "utopia_fe.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>


#include "libmesh/nonlinear_implicit_system.h"
#include "utopia_FEDSLMortarExamples.hpp"

#include "utopia_FEDSLMortarExamples.hpp"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/mesh_refinement.h"
#include "LibmeshTransferForMoose.hpp"
#include "LibmeshMultigridForMoose.hpp"


namespace libMesh {
	class LibMeshInit;
}

namespace utopia {

	using namespace std;
	using namespace libMesh;


  	void get_interpolation(libMesh::Parallel::Communicator &libmesh_comm, const std::shared_ptr<Mesh> &mesh_master, const std::shared_ptr<Mesh> &mesh_slave, DSMatrixd & T)
	{

        auto order_elem = FIRST;
        int order_quad = order_elem + order_elem;
        
        LibMeshFEContext<LinearImplicitSystem> master_context(mesh_master);
        auto master_space = fe_space(LAGRANGE, order_elem, master_context);
        master_context.equation_systems.init();
        
        LibMeshFEContext<LinearImplicitSystem> slave_context(mesh_slave);
        auto slave_space = fe_space(LAGRANGE, order_elem, slave_context);
        slave_context.equation_systems.init();
//        
//        MixedParMortarAssembler assembler(libmesh_comm, make_ref(master_space), make_ref(slave_space));
//        
//        DSMatrixd matrix;
//        assembler.Assemble_adaptive(matrix, 1, 0);
//        assembler.Transfer(matrix, T);
    }
    




	void linear_laplacian(LibMeshInit &init)
	{

		int n_fine  = 2;

		auto mesh_coarse = make_shared<Mesh>(init.comm());
		MeshTools::Generation::build_cube (*mesh_coarse,
			n_fine, n_fine, n_fine,
			-1., 1.,
			-1., 1.,
			-1., 1.,
			TET4);


		std::cout << "|master_elements| = " << mesh_coarse->n_elem() << std::endl;
		
		int coarse_elem = mesh_coarse->n_elem();


		EquationSystems master_es (*mesh_coarse);
		master_es.add_system<LinearImplicitSystem> ("TMP_MASTER");
		unsigned int u_var = master_es.get_system("TMP_MASTER").add_variable("v", FIRST);
		master_es.reinit();
		DofMap &master_dof = master_es.get_system(0).get_dof_map();

        
        std::vector<ElementDofMap> dof_map_coarse;
        std::vector<ElementDofMap> dof_map_fine;
        dof_map_coarse.resize(mesh_coarse->n_elem());
   
        dof_map_coarse.resize(mesh_coarse->n_elem());
        MeshBase::const_element_iterator e_it        =mesh_coarse->elements_begin();
        const MeshBase::const_element_iterator e_end =mesh_coarse->elements_end();
        std::vector<dof_id_type> temp;
        
    
        for (; e_it != e_end; ++e_it){
            
            Elem *elem = *e_it;
            master_dof.dof_indices(elem, temp, 0);
            
            dof_map_coarse[elem->id()].global.insert(dof_map_coarse[elem->id()].global.end(), temp.begin(), temp.end());
                
                std::cout<<"elem_id: "<< elem->id() << " tmp: "<< temp[0] << "  \n";
                std::cout<<"elem_id: "<< elem->id() << " tmp: "<< temp[1] << "  \n";
                std::cout<<"elem_id: "<< elem->id() << " tmp: "<< temp[2] << "  \n";
                std::cout<<"elem_id: "<< elem->id() << " tmp: "<< temp[3] << "  \n";

        }
        


     
      // auto mesh_coarse_copy= mesh_coarse;
//        
//        Mesh mesh_proof(init.comm());
//        
//        MeshTools::Generation::build_cube (mesh_proof,
//                                           n_fine, n_fine, n_fine,
//                                           -1., 1.,
//                                           -1., 1.,
//                                           -1., 1.,
//                                           TET4);
        
	    MeshRefinement mesh_refinement(*mesh_coarse);
      	mesh_refinement.uniformly_refine(1);


      	int fine_elem = mesh_coarse->n_elem() - coarse_elem;


      	std::cout << "|refined_elements| = " << mesh_coarse->n_elem() << std::endl;

		LibMeshFEContext<LinearImplicitSystem> context(mesh_coarse);
		auto Vh = fe_space(LAGRANGE, FIRST, context);
		auto v  = fe_function(Vh);

		strong_enforce( boundary_conditions(v == coeff(-1.0), {0, 1, 2, 3}) );

		context.equation_systems.init();

		const int dim = mesh_coarse->mesh_dimension();
		v.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));


		std::function<Real (const Point &)> f = [dim](const Point &p) {
			return 10 * std::sqrt( p(0) * p(0) + p(1) * p(1)) - 5.0;
		};

		auto ass = make_assembly([&]() -> void {
			std::cout << "assemble called, norm(u):" << context.system.solution->l2_norm() << std::endl;
			assemble(v, v, integral(dot(grad(v), grad(v)) ), integral(dot(coeff(f), v)), 
				*context.system.matrix,  
				*context.system.rhs);
		});



		context.system.attach_assemble_object(ass);
		context.equation_systems.print_info();
		// context.equation_systems.solve();

		// ExodusII_IO(*mesh_fine).write_equation_systems ("l_laplacian.e", context.equation_systems);




		auto ls_system   = dynamic_cast <LinearImplicitSystem *> (&context.equation_systems.get_system(0));  
		ls_system->assemble(); 

			long nn = ls_system->solution->size();
			DSMatrixd A;
			DVectord b; 
			DVectord sol; 
			sol = zeros(nn);
			A = sparse(nn, nn, 40);
		    b = zeros(nn);

	        Vec p_vec = cast_ptr< libMesh::PetscVector<libMesh::Number> *>(ls_system->rhs)->vec();
	        Mat p_mat = cast_ptr< libMesh::PetscMatrix<libMesh::Number> *>(ls_system->matrix)->mat();
	        VecAssemblyEnd(p_vec); 
	        MatAssemblyEnd(p_mat, MAT_FINAL_ASSEMBLY); 
			utopia::convert(p_vec, b);
			utopia::convert(p_mat, A);



	        int n_coarse  = 2;
//	       
//	        auto mesh_coarse = make_shared<Mesh>(init.comm());
//	        MeshTools::Generation::build_cube (*mesh_coarse,
//	                                           n_coarse, n_coarse, n_coarse,
//	                                           -1., 1.,
//	                                           -1., 1.,
//	                                           -1., 1.,
//	                                           TET4);
        
//           EquationSystems equation_systems(*mesh_coarse);
//           LinearImplicitSystem & proj_sys = equation_systems.add_system<LinearImplicitSystem>("TMP_slave");
//           unsigned int  _proj_var_num = proj_sys.add_variable("var", FIRST);
//           DofMap &coarse_dof = equation_systems.get_system("TMP_slave").get_dof_map();
        
	        using namespace utopia; 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



      	// slave
        context.equation_systems.add_system<LinearImplicitSystem> ("TMP_slave");
        context.equation_systems.get_system("TMP_slave").add_variable("u", FIRST);
        context.equation_systems.reinit(); 
        DofMap &slave_dof = context.equation_systems.get_system("TMP_slave").get_dof_map();

		// EquationSystems master_es (*mesh_coarse);
//		 master_es.add_system<LinearImplicitSystem> ("TMP_MASTER");
//		 unsigned int u_var = master_es.get_system("TMP_MASTER").add_variable("u", FIRST);
//		 master_es.reinit();

        
        

        
        dof_map_fine.resize(mesh_coarse->n_elem());
        MeshBase::const_element_iterator e_it_f        =mesh_coarse->local_level_elements_begin(1);
        const MeshBase::const_element_iterator e_end_f =mesh_coarse->local_level_elements_end(1);
        std::vector<dof_id_type> temp_f;
        
        
        for (; e_it_f != e_end_f; ++e_it_f){
            
            Elem *elem_f = *e_it_f;
            Elem *first_elem =*(mesh_coarse->local_level_elements_begin(1));
            slave_dof.dof_indices(elem_f, temp_f, 0);
            
            dof_map_fine[elem_f->id()-first_elem->id()].global.insert(dof_map_fine[elem_f->id()-first_elem->id()].global.end(), temp_f.begin(), temp_f.end());
            
            std::cout<<"elem_id: "<< elem_f->id()-first_elem->id() << " tmp: "<< temp_f[0] << "  \n";
            std::cout<<"elem_id: "<< elem_f->id()-first_elem->id() << " tmp: "<< temp_f[1] << "  \n";
            std::cout<<"elem_id: "<< elem_f->id()-first_elem->id() << " tmp: "<< temp_f[2] << "  \n";
            std::cout<<"elem_id: "<< elem_f->id()-first_elem->id() << " tmp: "<< temp_f[3] << "  \n";
            
        }

        
       // exit(1);
        
        

        unsigned int var_m=0, var_n=0;


		DSMatrixd _B; 


		// moonolith:Communicator expressComm(init.comm());
		moonolith:Communicator expressComm(PETSC_COMM_WORLD);



        // AssembleMOOSE_adaptive(expressComm, mesh_fine, 	mesh_fine,
        // 							utopia::make_ref(master_dof), 	utopia::make_ref(slave_dof), 
        // 							utopia::make_ref(var_m), 		utopia::make_ref(var_n), _B, 0, 1);



//		
//
//
//
		AssembleMultigridMOOSE(	expressComm,
                       			mesh_coarse,
		                       	utopia::make_ref(dof_map_coarse),
		                       	utopia::make_ref(dof_map_fine),
		                       	utopia::make_ref(var_m),
		                        utopia::make_ref(0),
		                       //const std::shared_ptr<const unsigned int> &master_elem, 
		                        utopia::make_ref(coarse_elem),
		                        utopia::make_ref(fine_elem),
		                       _B);







	    DVectord  diag_elem = 1./sum(_B,1);
	    DSMatrixd T = diag(diag_elem)*(_B);  		// interpolation operator for one variable


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//          DSMatrixd T; 
	//          get_interpolation(init.comm(), mesh_coarse, mesh_fine, T);

	//         get_interpolation(init.comm(), mesh_fine, mesh_fine, T);


	        std::cout<<"size:  "<< T.size().get(0)<< "   :   "<< T.size().get(1) << "  \n"; 


	        std::vector <DSMatrixd> interpolation_operators;
	        interpolation_operators.push_back(std::move(T));
	        
	        auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
	#ifdef PETSC_HAVE_MUMPS           
	        direct_solver->set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
	#endif //PETSC_HAVE_MUMPS            

	        auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
	        Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);


	        multigrid.init_transfer(std::move(interpolation_operators));
	        multigrid.galerkin_assembly(A);


	        Parameters params; 
	        params.verbose(true); 
	        params.linear_solver_verbose(true); 
	        multigrid.set_parameters(params); 

	        multigrid.solve(b, sol);







	}


	void run_solver_ex(libMesh::LibMeshInit &init)
	{
		linear_laplacian(init); 
	};
}

#endif //FEDSL_SOLVER_EXAMPLES_HPP
