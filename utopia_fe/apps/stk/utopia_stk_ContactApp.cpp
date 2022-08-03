
#include "utopia_Main.hpp"

#include "utopia_ContactApp.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

void stk_contact(utopia::Input &in) {
    utopia::ContactApp<utopia::stk::FunctionSpace> obs;
    obs.run(in);
}

UTOPIA_REGISTER_APP(stk_contact);
