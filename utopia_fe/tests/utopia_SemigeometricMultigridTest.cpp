#include "utopia_libmesh.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia.hpp"
#include "moonolith_communicator.hpp"
#include <memory>
#include "libmesh/parallel_mesh.h"
#include "libmesh/nemesis_io.h"

#include "utopia_SemiGeometricMultigrid.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_SemiGeometricMultigrid.hpp"

using namespace libMesh;

namespace utopia {

    void run_semigeometric_multigrid_elast(libMesh::Parallel::Communicator &comm)
    {
        std::cout << "[run_semigeometric_multigrid_elast]" << std::endl;

        auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(comm);

        const unsigned int n = 70;
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

        const double mu = 10;
        const double lambda = 10;

        auto e_u = 0.5 * ( transpose(grad(u)) + grad(u) );
        auto e_v = 0.5 * ( transpose(grad(v)) + grad(v) );

        // LMDenseVector z = zeros(2);
        LMDenseVector z = values(2, -0.2);
        auto elast_op = ((2. * mu) * inner(e_u, e_v) + lambda * inner(div(u), div(v))) * dX;
        auto f = inner(coeff(z), v) * dX;

        auto constr = constraints(
            boundary_conditions(uy == coeff(0.2),  {0}),
            // boundary_conditions(uy == coeff(-0.2),  {2}),
            boundary_conditions(uy == coeff(0.),  {2}),
            boundary_conditions(ux == coeff(0.0),  {0, 2})
            );

        init_constraints(constr);
        equation_systems->init();

        USparseMatrix stiffness_mat;
        UVector rhs;
        assemble(elast_op, stiffness_mat);
        assemble(f, rhs);
        apply_boundary_conditions(Vx.dof_map(), stiffness_mat, rhs);

        std::cout << "assembly complete" << std::endl;


        // auto linear_solver = std::make_shared<ConjugateGradient<USparseMatrix, UVector, HOMEMADE>>();
        // auto linear_solver = std::make_shared<BiCGStab<USparseMatrix, UVector>>();
        // auto linear_solver = std::make_shared<ConjugateGradient<USparseMatrix, UVector>>();
        auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
        auto smoother      = std::make_shared<GaussSeidel<USparseMatrix, UVector>>();
        // auto smoother = std::make_shared<ProjectedGaussSeidel<USparseMatrix, UVector>>();
        // auto smoother = std::make_shared<ConjugateGradient<USparseMatrix, UVector, HOMEMADE>>();
        // linear_solver->verbose(true);
        SemiGeometricMultigrid mg(smoother, linear_solver);
        // mg.algebraic().rtol(1e-16);
        // mg.algebraic().atol(1e-16);
        // mg.algebraic().max_it(400);

        // mg.convert_to_block_solver();
        mg.verbose(true);
        mg.init(sys, 4);
        // mg.max_it(1);

        UVector sol = local_zeros(local_size(rhs));

        Chrono c;
        c.start();
        mg.solve(stiffness_mat, rhs, sol);
        c.stop();

        if(mpi_world_rank() == 0) {
            std::cout << "multigrid solver:\n" << c << std::endl;
        }


        //CG with multigrid preconditioner
        // ConjugateGradient<USparseMatrix, UVector, HOMEMADE> cg;
        // BiCGStab<USparseMatrix, UVector> cg;
        // cg.verbose(true);
        // cg.set_preconditioner(make_ref(mg));
        // cg.solve(stiffness_mat, rhs, sol);

        const double err = norm2(stiffness_mat * sol - rhs);
        convert(sol, *sys.solution);
        sys.solution->close();
        Nemesis_IO(*lm_mesh).write_equation_systems("elast_mg.e", *equation_systems);
        assert(err < 1e-6);
    }

    void run_semigeometric_multigrid_poisson(libMesh::Parallel::Communicator &comm)
    {
        std::cout << "[run_semigeometric_multigrid_poisson]" << std::endl;
        auto lm_mesh = std::make_shared<libMesh::DistributedMesh>(comm);

        const unsigned int n = 50;
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

        USparseMatrix lapl_mat;
        UVector rhs;

        init_constraints(constr);
        equation_systems->init();

        assemble(lapl, lapl_mat);
        assemble(f, rhs);

        apply_boundary_conditions(V.dof_map(), lapl_mat, rhs);

         // auto linear_solver = std::make_shared<ConjugateGradient<USparseMatrix, UVector, HOMEMADE>>();
        auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
        // auto linear_solver = std::make_shared<ConjugateGradient<USparseMatrix, UVector>>();
        auto smoother = std::make_shared<GaussSeidel<USparseMatrix, UVector>>();
        // auto smoother = std::make_shared<ProjectedGaussSeidel<USparseMatrix, UVector, HOMEMADE>>();
        // smoother->set_n_local_sweeps(3);
        smoother->sweeps(3);

        // auto smoother = linear_solver;

        // linear_solver->verbose(true);
        SemiGeometricMultigrid mg(smoother, linear_solver);
        // mg.set_use_interpolation(true);
        mg.verbose(true);
        mg.init(V, 4);
        mg.update(make_ref(lapl_mat));

        UVector sol = local_zeros(local_size(rhs));

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

    void SMGTest::run(Input &in)
    {
        run_semigeometric_multigrid_poisson(comm());
        run_semigeometric_multigrid_elast(comm());
    }
}
