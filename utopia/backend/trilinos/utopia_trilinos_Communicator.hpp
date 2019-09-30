#ifndef UTOPIA_TRILINOS_COMMUNICATOR_HPP
#define UTOPIA_TRILINOS_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>

namespace utopia {

    class TrilinosCommunicator final : public Communicator {
    public:
        using CommPtr = Teuchos::RCP<const Teuchos::Comm<int> >;

        inline int rank() const override
        {
          return comm_->getRank();
        }

        inline int size() const override
        {
            return comm_->getSize();
        }

        inline Communicator * clone() const override
        {
            return new TrilinosCommunicator(get());
        }

        inline CommPtr get() const //override
        {
            return comm_;
        }

        inline void set(const CommPtr &comm)
        {
            comm_ = comm;
        }

        template<typename T>
        inline T sum(const T &val) const {

            T ret_global = 0.0;
            Teuchos::reduceAll(*get(), Teuchos::REDUCE_SUM, 1, &val, &ret_global);
            return ret_global;
        }

        TrilinosCommunicator(const CommPtr &comm) : comm_(comm) {}
        TrilinosCommunicator();

    private:
        CommPtr comm_;
    };
}


#endif //UTOPIA_TRILINOS_COMMUNICATOR_HPP
