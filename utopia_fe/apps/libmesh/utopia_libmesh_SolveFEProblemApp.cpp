
#include "utopia_Main.hpp"

#include "utopia_SolveFEProblemApp.hpp"

#include "utopia_libmesh.hpp"

void libmesh_solve(utopia::Input &in) {
    utopia::SolveFEProblemApp<utopia::libmesh::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(libmesh_solve);
