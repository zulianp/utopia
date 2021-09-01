#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_MortarApp.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::stk::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::stk::FunctionSpace>;
}  // namespace utopia

void stk_mortar(utopia::Input &in) {
    utopia::MortarApp<utopia::stk::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_mortar);
