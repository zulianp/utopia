#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_FracNetFlow.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"
#include "utopia_stk_intrepid2_Transport.hpp"

// REMOVE ME
#include "par_moonolith_instance.hpp"

void stk_fracflow(utopia::Input &in) {
    // REMOVE ME
    moonolith::Moonolith::instance().verbose(true);

    utopia::FracNetFlow<utopia::stk::FunctionSpace> fnf;
    fnf.read(in);  // Read from input file
    fnf.solve();   // Solve coupled problem
}

UTOPIA_REGISTER_APP(stk_fracflow);
