#include "utopia_trilinos_Communicator.hpp"
#include <Tpetra_Core.hpp>

namespace utopia {
    TrilinosCommunicator::TrilinosCommunicator()
    : comm_( Tpetra::getDefaultComm() )
    {

    }
}