#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_FracNetFlow.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_libmesh_kokkos_OmniAssembler.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"

// REMOVE ME
#include "par_moonolith_instance.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::libmesh::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::libmesh::FunctionSpace>;
}  // namespace utopia

void libmesh_fracflow(utopia::Input &in) {
    // REMOVE ME
    moonolith::Moonolith::instance().verbose(true);

    utopia::FracNetFlow<utopia::libmesh::FunctionSpace> fnf;
    fnf.read(in);
    fnf.solve();
}

UTOPIA_REGISTER_APP(libmesh_fracflow);
