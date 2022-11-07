#include "utopia_Main.hpp"

#include "utopia.hpp"
#include "utopia_kokkos_FunctionSpace.hpp"

// #include "utopia_SimpleNewton.hpp"

// #include "utopia_ImplicitEulerIntegrator.hpp"
// #include "utopia_NLSolveApp.hpp"
// #include "utopia_NewmarkIntegrator.hpp"
// #include "utopia_SemiGeometricMultigridNew.hpp"

// #include "utopia_kokkos.hpp"
// #include "utopia_kokkos_OmniAssembler.hpp"
// // #include "utopia_moonolith_kokkos_FETransfer.hpp"

// namespace utopia {
//     template class NewmarkIntegrator<utopia::kokkos::FunctionSpace>;
//     template class ImplicitEulerIntegrator<utopia::kokkos::FunctionSpace>;
// }  // namespace utopia

void kokkos_nlsolve(utopia::Input &in) {
    // utopia::NLSolveApp<utopia::kokkos::FunctionSpace> app;
    // app.read(in);

    // if (app.valid()) {
    //     app.run();
    // } else {
    //     utopia::err() << "kokkos_nlsolve: invalid app setup\n";
    // }

    using Scalar = double;
    utopia::kokkos::FunctionSpace<Scalar> space;
    in.get("space", space);
}

UTOPIA_REGISTER_APP(kokkos_nlsolve);
