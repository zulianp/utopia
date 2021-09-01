#include "utopia_Main.hpp"

#include "utopia_FETransferApp.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_FSIHesch2014.hpp"

void stk_fsi(utopia::Input &in) {
    using FunctionSpace = utopia::stk::FunctionSpace;
    using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;

    std::shared_ptr<FEFunctionInterface_t> fluid, solid;
    utopia::FSIHesch2014<FunctionSpace> fsi(fluid, solid);
    fsi.read(in);
}

UTOPIA_REGISTER_APP(stk_fsi);
