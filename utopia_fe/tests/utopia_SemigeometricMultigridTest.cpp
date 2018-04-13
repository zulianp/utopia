#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia.hpp"
#include "moonolith_communicator.hpp"
#include <memory>
#include "libmesh/parallel_mesh.h"
#include "libmesh/nemesis_io.h"

#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

using namespace libMesh;

namespace utopia {

    void init_mg_test_problem(libMesh::LibMeshInit &init, MGTestProblem &p)
    {
        using std::make_shared;
        
        Order order_elem_fine     = FIRST;
        Order order_elem_coarse   = FIRST;
        Order order_of_quadrature = Order(int(order_elem_fine) + int(order_elem_coarse));
        
        const std::string data_path = Utopia::instance().get("data_path");
        const std::string fine_mesh_path   = data_path + "/fine_mesh.e";
        const std::string coarse_mesh_path = data_path + "/coarse_mesh.e";
        
        p.fine_mesh = make_shared<DistributedMesh>(init.comm());
        p.fine_mesh->read(fine_mesh_path);
        
        p.coarse_mesh = make_shared<DistributedMesh>(init.comm());
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
        
//        if(!assemble_volume_transfer(comm,
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
//        if(!assemble_volume_transfer(comm,
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
        if(!assemble_volume_transfer_r(comm,
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
         B_r))
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
        // p.fine_context->equation_systems.print_info();
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
    }
    
    void run_semigeometric_multigrid_test_old(libMesh::LibMeshInit &init)
    {
        MGTestProblem p;
        init_mg_test_problem(init, p);
        DVectord sol = local_zeros(local_size(p.rhs));
        
        //set up multigrid
        // auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord> >();
        auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
        auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
        
        Multigrid<DSMatrixd, DVectord> multigrid(smoother, linear_solver);
        multigrid.init_transfer_from_coarse_to_fine(p.interpolation_operators);
        multigrid.mg_type(2);
        
        // multigrid.max_it(1);
        // multigrid.verbose(true);
        multigrid.solve(p.A, p.rhs, sol);

        const double err = norm2(p.A * sol - p.rhs);
        
        
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
        assert(err < 1e-7);
    }

    //not working properly
    void run_semigeometric_multigrid_elast(libMesh::LibMeshInit &init)
    {
        std::cout << "[run_semigeometric_multigrid_elast]" << std::endl;

        auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());     
        
        const unsigned int n = 100;
        libMesh::MeshTools::Generation::build_square(*lm_mesh,
            n, n,
            0, 1,
            0, 1.,
            libMesh::QUAD4);

        int dim = lm_mesh->mesh_dimension();

        auto equation_systems = std::make_shared<libMesh::EquationSystems>(*lm_mesh);
        auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("smg_elast");

        auto Vx = LibMeshFunctionSpace(equation_systems);
        auto Vy = LibMeshFunctionSpace(equation_systems);
        auto V = Vx * Vy;

        auto u = trial(V);
        auto v = test(V);

        auto ux = u[0];
        auto uy = u[1];

        const double mu = 1;
        const double lambda = 1;

        auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) ); 
        auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

        // LMDenseVector z = zeros(2);
        LMDenseVector z = values(2, -0.2);
        auto elast_op = ((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v))) * dX;
        auto f = inner(coeff(z), v) * dX;

        auto constr = constraints(
            boundary_conditions(uy == coeff(0.2),  {0}),
            boundary_conditions(uy == coeff(0.0),  {2}),
            boundary_conditions(ux == coeff(0.0),  {0, 2})
            );

        init_constraints(constr);
        equation_systems->init();

        DSMatrixd stiffness_mat;
        DVectord rhs;
        assemble(elast_op, stiffness_mat);
        assemble(f, rhs);     
        apply_boundary_conditions(Vx.dof_map(), stiffness_mat, rhs);

        std::cout << "assembly complete" << std::endl;

        auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>();
        // auto linear_solver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>();
        // auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
        // auto linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>();
        // auto smoother      = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
        // auto smoother = std::make_shared<ProjectedGaussSeidel<DSMatrixd, DVectord>>();
        auto smoother = std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>();
        // linear_solver->verbose(true);
        SemiGeometricMultigrid mg(smoother, linear_solver);
        mg.algebraic().rtol(1e-16);
        mg.algebraic().atol(1e-16);
        mg.algebraic().max_it(400);

        // mg.convert_to_block_solver();
        mg.verbose(true);
        mg.init(*equation_systems, 4);
        // mg.max_it(1);

        DVectord sol = local_zeros(local_size(rhs));

        Chrono c;
        c.start();
        mg.solve(stiffness_mat, rhs, sol);
        c.stop();
        
        if(mpi_world_rank() == 0) {
            std::cout << "multigrid solver:\n" << c << std::endl;
        }


        //CG with multigrid preconditioner
        // ConjugateGradient<DSMatrixd, DVectord, HOMEMADE> cg;
        // BiCGStab<DSMatrixd, DVectord> cg;
        // cg.verbose(true);
        // cg.set_preconditioner(make_ref(mg));
        // cg.solve(stiffness_mat, rhs, sol);

        const double err = norm2(stiffness_mat * sol - rhs);
        convert(sol, *sys.solution);
        sys.solution->close();
        Nemesis_IO(*lm_mesh).write_equation_systems("elast_mg.e", *equation_systems);
        assert(err < 1e-6);
    }

    void run_semigeometric_multigrid_poisson(libMesh::LibMeshInit &init)
    {
        std::cout << "[run_semigeometric_multigrid_poisson]" << std::endl;
        auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());     
        
        const unsigned int n = 200;
        libMesh::MeshTools::Generation::build_square(*lm_mesh,
            n, n,
            0, 1,
            0, 1.,
            libMesh::QUAD4);

        auto equation_systems = std::make_shared<libMesh::EquationSystems>(*lm_mesh);
        equation_systems->add_system<libMesh::LinearImplicitSystem>("smg");

        auto V = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, libMesh::FIRST, "u");
        auto u = trial(V);
        auto v = test(V);

        auto lapl = inner(grad(u), grad(v)) * dX;
        auto f = inner(coeff(1.), v) * dX;

        auto constr = constraints(
            boundary_conditions(u == coeff(0.),  {1, 3}),
            boundary_conditions(u == coeff(0.),  {0}),
            boundary_conditions(u == coeff(0.0), {2})
            );

        DSMatrixd lapl_mat;
        DVectord rhs;

        init_constraints(constr);
        equation_systems->init();

        assemble(lapl, lapl_mat);
        assemble(f, rhs);

        apply_boundary_conditions(V.dof_map(), lapl_mat, rhs);

         // auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord, HOMEMADE>>();
        // auto linear_solver = std::make_shared<BiCGStab<DSMatrixd, DVectord>>());
        auto linear_solver = std::make_shared<ConjugateGradient<DSMatrixd, DVectord>>();
        // auto smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>();
        auto smoother = std::make_shared<ProjectedGaussSeidel<DSMatrixd, DVectord, HOMEMADE>>();
        smoother->set_n_local_sweeps(3);
        smoother->sweeps(1);
        
        // auto smoother = linear_solver;

        // linear_solver->verbose(true);
        SemiGeometricMultigrid mg(smoother, linear_solver);
        mg.verbose(true);
        mg.init(V, 4);
        mg.update(make_ref(lapl_mat));

        DVectord sol = local_zeros(local_size(rhs));

        Chrono c;
        c.start();
        mg.apply(rhs, sol);
        c.stop();

        if(mpi_world_rank() == 0) {
            std::cout << "multigrid solver:\n" << c << std::endl;
        }

        const double err = norm2(lapl_mat * sol - rhs);
        assert(err < 1e-6);
    }

    void run_semigeometric_multigrid_test(libMesh::LibMeshInit &init)
    {
        run_semigeometric_multigrid_poisson(init);
        run_semigeometric_multigrid_elast(init);
    }
}
