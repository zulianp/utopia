#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NLSolveApp.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_mars.hpp"
#include "utopia_mars_OmniAssembler.hpp"
// #include "utopia_moonolith_mars_FETransfer.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::mars::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::mars::FunctionSpace>;
}  // namespace utopia

void mars_nlsolve(utopia::Input &in) {
    utopia::NLSolveApp<utopia::mars::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "mars_nlsolve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(mars_nlsolve);
