
#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_INTREPID2
#ifdef UTOPIA_WITH_PETSC_DM

#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NLSolveApp.hpp"
#include "utopia_NewmarkIntegrator.hpp"
// #include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_moonolith_petsc_FETransfer.hpp"
#include "utopia_petsc_dm.hpp"


#include "utopia_petsc_intrepid2_OmniAssembler.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::petsc::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::petsc::FunctionSpace>;
}  // namespace utopia

void petsc_nlsolve(utopia::Input &in) {
    utopia::NLSolveApp<utopia::petsc::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "petsc_nlsolve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(petsc_nlsolve);

#endif  //UTOPIA_WITH_INTREPID2
#endif  //UTOPIA_WITH_PETSC_DM