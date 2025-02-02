#include "utopia_fe_base.hpp"

#ifdef UTOPIA_ENABLE_INTREPID2
#ifdef UTOPIA_ENABLE_PETSC_DM

#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_FraNetGApp.hpp"

#include "utopia_moonolith_petsc_FETransfer.hpp"
#include "utopia_petsc_dm.hpp"

#include "utopia_petsc_intrepid2_OmniAssembler.hpp"

void petsc_franetg(utopia::Input &in) {
    utopia::FraNetGApp<utopia::petsc::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(petsc_franetg);

#endif  // UTOPIA_ENABLE_INTREPID2
#endif  // UTOPIA_ENABLE_PETSC_DM