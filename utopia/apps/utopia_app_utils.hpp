#ifndef UTOPIA_APP_UTILS_HPP
#define UTOPIA_APP_UTILS_HPP

#include "utopia_GeometricMultigrid.hpp"
#include "utopia_MPITimeStatistics.hpp"

// #include "petscdraw.h"

namespace utopia {

    template<class Model, class FunctionSpace>
    static void geometric_multigrid(
        FunctionSpace &coarse_space,
        Input &in)
    {
        using Matrix = typename FunctionSpace::Matrix;
        using Vector = typename FunctionSpace::Vector;

        auto &comm = coarse_space.comm();
        MPITimeStatistics stats(comm); stats.start();

        int n_levels = 2;
        in.get("n_levels", n_levels);

        // bool in_situ_rendering = false;
        // in.get("in_situ_rendering", in_situ_rendering);

        std::string output_path = "MG.vtr";
        in.get("output_path", output_path);

        bool export_system = false;
        in.get("export_system", export_system);

        bool direct_solution = false;
        in.get("direct_solution", direct_solution);

    	Vector x;

        std::shared_ptr<FunctionSpace> fine_space_ptr;

    	{
        	auto smoother      = std::make_shared<SOR<Matrix, Vector>>();
            std::shared_ptr< LinearSolver<Matrix, Vector> > coarse_solver;
            // if(direct_solution) {
                coarse_solver = std::make_shared<Factorization<Matrix, Vector>>();
		    // } else {
      //           auto bcg = std::make_shared<BiCGStab<Matrix, Vector>>("bjacobi");
      //           // auto bcg = std::make_shared<GMRES<Matrix, Vector>>("bjacobi");
      //           bcg->max_it(coarse_space.n_dofs());
      //           coarse_solver = bcg;
      //       }

        	GeometricMultigrid<FunctionSpace> mg(smoother, coarse_solver);
        	// mg.verbose(true);
            mg.read(in);
        	mg.init(coarse_space, n_levels);

       		UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after mg-setup");
        	stats.stop_collect_and_restart("mg-setup");

        	Vector b;
        	Matrix A;

            fine_space_ptr = mg.fine_space_ptr();
        	auto &space = *fine_space_ptr;

        	space.create_vector(x);
        	space.create_vector(b);
        	space.create_matrix(A);

        	UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-create-matrix");

        	Model model(space);
        	model.read(in);

        	model.hessian(x, A);

        	UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-assembly");
        	stats.stop_collect_and_restart("tensor-creation+assembly");

        	space.apply_constraints(b);

            if(export_system) {
                rename("a", A);
                write("A.m", A);

                rename("b", b);
                write("B.m", b);
            }

            if(direct_solution) {
                Factorization<Matrix, Vector> solver;
                solver.solve(A, b, x);
                rename("x", x);
                fine_space_ptr->write("direct_solution.vtr", x);
                x.set(0.0);
            }

            ConjugateGradient<Matrix, Vector, HOMEMADE> cg;
            // BiCGStab<Matrix, Vector, HOMEMADE> cg;
            cg.read(in);

            mg.max_it(1);
            cg.set_preconditioner(make_ref(mg));

            cg.update(make_ref(A));
        	UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-update");

            cg.verbose(true);

            if(!cg.apply(b, x)) {
                std::cerr << "[Error] Unable to solve system up to requested precision" << std::endl;
            }

        	UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-apply");
        	stats.stop_collect_and_restart("solve");
    	}

        if(!output_path.empty()) {
            rename("x", x);
            fine_space_ptr->write(output_path, x);
            stats.stop_collect_and_restart("output");
    	}

        comm.root_print("n_dofs: " + std::to_string(x.size()));
        stats.describe(std::cout);
    }

}


#endif //UTOPIA_APP_UTILS_HPP
