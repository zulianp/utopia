#ifndef UTOPIA_PETSC_COMMUNICATOR_HPP
#define UTOPIA_PETSC_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Utils.hpp"

#include <mpi.h>
#include <memory>

namespace utopia {

    class PetscCommunicator final : public MPICommunicator {
    public:
        inline PetscCommunicator *clone() const override { return new PetscCommunicator(get()); }

        static PetscCommunicator self();
        static PetscCommunicator world();

        inline static PetscCommunicator &get_default() {
            static PetscCommunicator instance_;
            return instance_;
        }

        PetscCommunicator();

        explicit PetscCommunicator(const MPI_Comm comm) : MPICommunicator(comm, false) {}
        explicit PetscCommunicator(const Communicator &comm) : MPICommunicator(comm) {}

        inline PetscCommunicator(const MPICommunicator &other) : MPICommunicator(other) {}
        inline PetscCommunicator(MPICommunicator &&other) : MPICommunicator(std::move(other)) {}
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_COMMUNICATOR_HPP
