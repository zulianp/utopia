#include "utopia_Main.hpp"

#include "utopia_FETransferTimeSeriesApp.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"
#include "utopia_stk_intrepid2.hpp"

void stk_transfer_time_series(utopia::Input &in) {
    utopia::FETransferTimeSeriesApp<utopia::stk::FunctionSpace> transfer;
    transfer.read(in);
    if (transfer.is_valid()) {
        transfer.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(stk_transfer_time_series);
