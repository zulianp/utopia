#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS

// include edsl components
#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_petsc_trilinos.hpp"
#include "utopia_trilinos.hpp"

#include "utopia_MPITimeStatistics.hpp"

// std
#include <cmath>

namespace utopia {

    void ls_solve(Input &in) {
        TpetraMatrix A;
        TpetraVector x, b, sol;

        MPITimeStatistics stats(x.comm());

        //////////////////////////////////////////
        stats.start();

        read("../data/knf/matrices/NavierStokes-O-0.mm.1", A);
        read("../data/knf/matrices/NavierStokes-O-0.rhs.1", b);
        read("../data/knf/matrices/NavierStokes-O-0.sln.1", sol);

        stats.stop_collect_and_restart("read_files");

        // disp(b);

        // BelosSolver<TpetraMatrix, TpetraVector> solver;
        KSPSolver<TpetraMatrix, TpetraVector> solver;
        solver.read(in);

        stats.stop_collect_and_restart("read_settings");

        disp("solve begin");
        solver.solve(A, b, x);
        disp("solve end");

        stats.stop_collect_and_restart("solve");

        write("../data/knf/matrices/result.mm", x);

        stats.stop_collect_and_restart("write");
        stats.describe(utopia::out().stream());
    }

    UTOPIA_REGISTER_APP(ls_solve);

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS
