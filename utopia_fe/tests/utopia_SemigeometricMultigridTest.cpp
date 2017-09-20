#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia.hpp"
#include "moonolith_communicator.hpp"
#include <memory>

using namespace libMesh;

namespace utopia {
    
    void init_mg_test_problem(libMesh::LibMeshInit &init, MGTestProblem &p)
    {
        using std::make_shared;
        
        Order order_elem_fine     = FIRST;
        Order order_elem_coarse   = FIRST;
        Order order_of_quadrature = Order(int(order_elem_fine) + int(order_elem_coarse));
        
        const std::string data_path = Utopia::Instance().get("data_path");
        const std::string fine_mesh_path   = data_path + "/fine_mesh.e";
        const std::string coarse_mesh_path = data_path + "/coarse_mesh.e";
        
        p.fine_mesh = make_shared<Mesh>(init.comm());
        p.fine_mesh->read(fine_mesh_path);
        
        p.coarse_mesh = make_shared<Mesh>(init.comm());
        p.coarse_mesh->read(coarse_mesh_path);
        
        p.fine_context 	 = make_shared<MGTestProblem::FEContextT>(p.fine_mesh);
        p.coarse_context = make_shared<MGTestProblem::FEContextT>(p.coarse_mesh);
        
        p.fine_space   = make_shared<MGTestProblem::FESpaceT>( fe_space(LAGRANGE, order_elem_fine, *p.fine_context) );
        p.coarse_space = make_shared<MGTestProblem::FESpaceT>( fe_space(LAGRANGE, order_elem_coarse, *p.coarse_context) );
        
        auto u = fe_function(*p.fine_space);
        strong_enforce( boundary_conditions(u == coeff(0.), {1}) );
        
        //this needs to be initialized here
        p.fine_context->equation_systems.init();
        p.coarse_context->equation_systems.init();
        
        moonolith::Communicator comm(init.comm().get());
        DSMatrixd B;
                  /* without volume tags*/
        
//        if(!assemble_volume_transfer(
//                                     comm,
//                                     p.coarse_mesh,
//                                     p.fine_mesh,
//                                     utopia::make_ref(p.coarse_space->dof_map()),
//                                     utopia::make_ref(p.fine_space->dof_map()),
//                                     0,
//                                     0,
//                                     true,
//                                     1,
//                                     B))
        
                  /* with volume tags*/
            
//        if(!assemble_volume_transfer(
//                                     comm,
//                                     p.coarse_mesh,
//                                     p.fine_mesh,
//                                     utopia::make_ref(p.coarse_space->dof_map()),
//                                     utopia::make_ref(p.fine_space->dof_map()),
//                                     0,
//                                     0,
//                                     true,
//                                     1,
//                                     B,
//                                     {{1,2}}))
        
                  /* with volume tags and reverse operator*/
        DSMatrixd B_r;
        if(!assemble_volume_transfer_r(
                                       comm,
                                       p.coarse_mesh,
                                       p.fine_mesh,
                                       utopia::make_ref(p.coarse_space->dof_map()),
                                       utopia::make_ref(p.fine_space->dof_map()),
                                       utopia::make_ref(p.coarse_space->dof_map()),
                                       utopia::make_ref(p.fine_space->dof_map()),
                                       0,
                                       0,
                                       0,
                                       0,
                                       true,
                                       1,
                                       1,
                                       B,
                                       B_r
                                       ))
        {
            
            std::cerr << "No intersection" << std::endl;
            return;
        }
        
        const int dim = p.fine_mesh->mesh_dimension();
        
        auto b_form = integral(dot(grad(u), grad(u)));
        auto l_form = integral(dot(coeff(1.0), u));
        
        u.set_quad_rule(std::make_shared<libMesh::QGauss>(dim, order_of_quadrature));
        
        auto ass = make_assembly([&]() -> void {
            assemble(u, u, b_form, l_form, *p.fine_context->system.matrix, *p.fine_context->system.rhs);
        });
        
        p.fine_context->system.attach_assemble_object(ass);
        p.fine_context->equation_systems.print_info();
        // p.fine_context->equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
        //solve triggers the initialization and the assembly of the algebraic system (should be fixed at some point)
        p.fine_context->equation_systems.solve();
        
        DVectord d = sum(B, 1);
        DVectord d_inv = local_values(local_size(d).get(0), 1.);
        
        {
            Write<DVectord> w(d_inv);
            each_read(d, [&](const SizeType i, const double value) {
                if(value > 1e-15) {
                    d_inv.set(i, 1./value);
                }
            });
        }
        
        DSMatrixd D_inv = diag(d_inv);
        
        std::shared_ptr<DSMatrixd> T = std::make_shared<DSMatrixd>();
        *T = D_inv * B;
        
        p.interpolation_operators.push_back(T);
        p.rhs = local_zeros(local_size(d));
        
        convert(*p.fine_context->system.matrix, p.A);
        convert(*p.fine_context->system.rhs, 	p.rhs);
        
        // write("T.m", *T);
    }
    
    void run_semigeometric_multigrid_test(libMesh::LibMeshInit &init)
    {
        MGTestProblem p;
        init_mg_test_problem(init, p);
        DVectord sol = local_zeros(local_size(p.rhs));
        
        //set up multigrid
        auto direct_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
        auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
        
        Multigrid<DSMatrixd, DVectord> multigrid(smoother, direct_solver);
        multigrid.init_transfer_from_coarse_to_fine(p.interpolation_operators);
        multigrid.mg_type(2);
        
        // multigrid.max_it(1);
        multigrid.verbose(true);
        multigrid.solve(p.A, p.rhs, sol);
        
        //CG with multigrid preconditioner
        // ConjugateGradient<DSMatrixd, DVectord> cg;
        // cg.verbose(true);
        // cg.set_preconditioner(make_ref(multigrid));
        // cg.solve(p.A, p.rhs, sol);
        
        DVectord ref_sol = local_zeros(local_size(sol));
        convert(*p.fine_context->system.solution, ref_sol);
        assert(approxeq(ref_sol, sol));
        
        convert(sol, *p.fine_context->system.solution);
        ExodusII_IO(*p.fine_mesh).write_equation_systems ("mg_solution_lapl.e", p.fine_context->equation_systems);
        
        
    }
}
