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
    using Scalar = double;

    utopia::MPICommunicator comm(MPI_COMM_WORLD);
    utopia::kokkos::FunctionSpace<Scalar> space(comm);
    in.get("space", space);
}

UTOPIA_REGISTER_APP(kokkos_nlsolve);
