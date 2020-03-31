#include "utopia_trilinos_Communicator.hpp"
#include <Tpetra_Core.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

namespace utopia {
    TrilinosCommunicator::TrilinosCommunicator()
    : comm_( Tpetra::getDefaultComm() )
    {

    }

   // TrilinosCommunicator::TrilinosCommunicator(const SelfCommunicator &comm) : comm_(new Teuchos::MpiComm<int>(comm.get())) {}

    TrilinosCommunicator::TrilinosCommunicator(const SelfCommunicator &comm) : comm_(new Teuchos::SerialComm<int>()) {}
}