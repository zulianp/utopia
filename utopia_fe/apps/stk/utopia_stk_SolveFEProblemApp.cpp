
#include "utopia_Main.hpp"

#include "utopia_SolveFEProblemApp.hpp"

#include "utopia_stk_intrepid2_OmniAssembler.hpp"

void stk_solve(utopia::Input &in) {
    utopia::SolveFEProblemApp<utopia::stk::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_solve);
