#include "utopia_Main.hpp"

#include "utopia_intrepid2_MechPostProcessApp.hpp"

#include "utopia_stk.hpp"

// Include stk to intrepid2 adapters
#include "utopia_stk_intrepid2.hpp"

void stk_mech_post_process(utopia::Input &in) {
    utopia::intrepid2::MechPostProcessApp<utopia::stk::FunctionSpace> grad;
    grad.read(in);
    if (grad.valid()) {
        grad.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(stk_mech_post_process);
