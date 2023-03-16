#include "utopia_BDDSQPApp.hpp"

#include "utopia_Main.hpp"

#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_moonolith_stk_Contact.hpp"
#include "utopia_moonolith_stk_FETransfer.hpp"
// #include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_stk_intrepid2_Discretization.hpp"
#include "utopia_stk_intrepid2_Material.hpp"

void stk_bddsqp(utopia::Input &in) {
    utopia::BDDSQPApp<utopia::stk::FunctionSpace, utopia::intrepid2::FE<double>> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_bddsqp);
