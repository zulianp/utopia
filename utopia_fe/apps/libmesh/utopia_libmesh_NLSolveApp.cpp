#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NLSolveApp.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_OmniAssembler.hpp"

#include "utopia_moonolith_libmesh_Contact.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"
#include "utopia_moonolith_libmesh_Obstacle.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::libmesh::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::libmesh::FunctionSpace>;
}  // namespace utopia

void libmesh_nlsolve(utopia::Input &in) {
    utopia::NLSolveApp<utopia::libmesh::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "libmesh_nlsolve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(libmesh_nlsolve);
