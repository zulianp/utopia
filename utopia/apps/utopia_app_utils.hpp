#ifndef UTOPIA_APP_UTILS_HPP
#define UTOPIA_APP_UTILS_HPP

#include "utopia_GeometricMultigrid.hpp"
#include "utopia_MPITimeStatistics.hpp"

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

        std::string output_path = "MG.vtr";
        in.get("output_path", output_path);


        auto smoother      = std::make_shared<SOR<Matrix, Vector>>();
        auto coarse_solver = std::make_shared<Factorization<Matrix, Vector>>();

        GeometricMultigrid<FunctionSpace> mg(smoother, coarse_solver);
        mg.verbose(true);
        mg.init(coarse_space, n_levels);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after mg-setup");
        stats.stop_collect_and_restart("mg-setup");

        Vector x, b;
        Matrix A;

        auto &space = mg.fine_space();

        space.create_vector(x);
        space.create_vector(b);
        space.create_matrix(A);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-create-matrix");

        Model model(space);
        model.read(in);

        model.hessian(x, A);

        // bool write_mat = false;
        // in.get("write_mat", write_mat);

        // if(write_mat) {
        //     rename("a", A);
        //     write("A.m", A);
        // }
        // model.gradient(x, b);

        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-assembly");
        stats.stop_collect_and_restart("tensor-creation+assembly");

        space.apply_constraints(b);

        mg.update(make_ref(A));
        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-update");

        mg.apply(b, x);
        UTOPIA_PETSC_COLLECTIVE_MEMUSAGE("after-apply");
        stats.stop_collect_and_restart("solve");

        rename("x", x);
        space.write(output_path, x);

        stats.stop_collect_and_restart("output");

        comm.root_print("n_dofs: " + std::to_string(space.n_dofs()));
        stats.describe(std::cout);
    }

}


#endif //UTOPIA_APP_UTILS_HPP