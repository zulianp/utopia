#ifndef UTOPIA_PETSC_COMMUNICATOR_HPP
#define UTOPIA_PETSC_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Utils.hpp"

#include <mpi.h>
#include <memory>

namespace utopia {

    class PetscCommunicator final : public MPICommunicator {
    public:
        class Wrapper {
        public:
            Wrapper(const MPI_Comm &comm, const bool owned) : comm(comm), owned(owned) {}

            ~Wrapper() {
                if (owned) {
                    MPI_Comm_free(&comm);
                }
            }

            MPI_Comm comm;
            bool owned;
        };

        inline PetscCommunicator *clone() const override { return new PetscCommunicator(get()); }

        inline MPI_Comm get() const override { return wrapper_->comm; }

        inline MPI_Comm raw_comm() const override { return get(); }

        inline void set(MPI_Comm comm) { wrapper_ = make_not_owned(comm); }

        inline void own(MPI_Comm comm) { wrapper_ = std::make_shared<Wrapper>(comm, true); }

        static PetscCommunicator self();
        static PetscCommunicator world();

        inline static PetscCommunicator &get_default() {
            static PetscCommunicator instance_;
            return instance_;
        }

        PetscCommunicator split(const int color) const;

        explicit PetscCommunicator(const MPI_Comm comm) : wrapper_(make_not_owned(comm)) {}
        PetscCommunicator();

        PetscCommunicator(const Communicator &comm) : wrapper_(make_not_owned(comm.raw_comm())) {}

        PetscCommunicator(const PetscCommunicator &other) : wrapper_(other.wrapper_) {}

    private:
        std::shared_ptr<Wrapper> wrapper_;

        static std::shared_ptr<Wrapper> make_not_owned(MPI_Comm comm) { return std::make_shared<Wrapper>(comm, false); }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_COMMUNICATOR_HPP
