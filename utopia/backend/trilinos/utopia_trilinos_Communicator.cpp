#include "utopia_Base.hpp"

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Tpetra_Core.hpp>
#include "utopia_trilinos_Communicator.hpp"

namespace utopia {
    TrilinosCommunicator::TrilinosCommunicator() : comm_(Tpetra::getDefaultComm()) {}

    // TrilinosCommunicator::TrilinosCommunicator(const SelfCommunicator &comm) : comm_(new
    // Teuchos::MpiComm<int>(comm.get())) {}

    TrilinosCommunicator::TrilinosCommunicator(const SelfCommunicator & /*comm*/)
        : comm_(new Teuchos::SerialComm<int>()) {}

    TrilinosCommunicator::TrilinosCommunicator(const Communicator &comm)
        : comm_(new Teuchos::MpiComm<int>(comm.raw_comm())) {
        assert(comm.size() == this->size());
        assert(comm.rank() == this->rank());
    }

    TrilinosCommunicator TrilinosCommunicator::split(const int color) const {
        return TrilinosCommunicator(get()->split(color, rank()));
    }

#ifdef UTOPIA_WITH_MPI
    MPI_Comm TrilinosCommunicator::raw_comm() const {
        auto *mpi_comm = dynamic_cast<const Teuchos::MpiComm<int> *>(comm_.get());
        if (mpi_comm != nullptr) {
            return *mpi_comm->getRawMpiComm();
        }
        assert(false);
        return MPI_COMM_SELF;
    }

    TrilinosCommunicator::TrilinosCommunicator(MPI_Comm comm) : comm_(new Teuchos::MpiComm<int>(comm)) {}

    TrilinosCommunicator TrilinosCommunicator::self() {
        return TrilinosCommunicator(Teuchos::rcp(new Teuchos::SerialComm<int>()));
    }

    TrilinosCommunicator TrilinosCommunicator::world() {
        // FIXME
        return TrilinosCommunicator(Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD)));
    }
#endif
}  // namespace utopia