#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NLSolveApp.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"

#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"

#ifdef UTOPIA_ENABLE_SFEM
#include "utopia_sfem_stk_SDFObstacle.hpp"
#endif  // UTOPIA_ENABLE_SFEM

void stk_nlsolve(utopia::Input &in) {
#ifdef UTOPIA_ENABLE_SFEM
    utopia::register_sfem_stk_contact();
#endif  // UTOPIA_ENABLE_SFEM

    utopia::NLSolveApp<utopia::stk::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "stk_nlsolve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(stk_nlsolve);
