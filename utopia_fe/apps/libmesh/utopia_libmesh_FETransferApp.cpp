
#include "utopia_Main.hpp"

#include "utopia_FETransferApp.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_moonolith_libmesh_FETransfer.hpp"

void libmesh_transfer(utopia::Input &in) {
    utopia::FETransferApp<utopia::libmesh::FunctionSpace> transfer;
    transfer.read(in);
    if (transfer.is_valid()) {
        transfer.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(libmesh_transfer);
