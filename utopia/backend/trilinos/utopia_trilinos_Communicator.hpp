#ifndef UTOPIA_TRILINOS_COMMUNICATOR_HPP
#define UTOPIA_TRILINOS_COMMUNICATOR_HPP

#include "utopia_Communicator.hpp"
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>

#include <array>

namespace utopia {

    class TrilinosCommunicator final : public Communicator {
    public:
        using CommPtr = Teuchos::RCP<const Teuchos::Comm<int> >;

        inline void barrier() const override
        {
            comm_->barrier();
        }

        inline int rank() const override
        {
          return comm_->getRank();
        }

        inline int size() const override
        {
            return comm_->getSize();
        }

#ifdef WITH_MPI
        MPI_Comm raw_comm() const override;
#endif

        inline TrilinosCommunicator * clone() const override
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

        inline bool conjunction(const bool &val) const override  {
            int i = val;
            i = sum(i);
            return i == this->size();
        }

        inline bool disjunction(const bool &val) const override  {
            int i = val;
            i = sum(i);
            return i > 0;
        }

        template<typename T>
        inline T sum(const T &val) const {

            T ret_global = 0.0;
            Teuchos::reduceAll(*get(), Teuchos::REDUCE_SUM, 1, &val, &ret_global);
            return ret_global;
        }

        template<typename T, std::size_t N>
        inline T sum(std::array<T, N> &val) const {
            T ret_global = 0.0;

            std::array<T, N> local_val;
            std::copy(std::begin(val), std::end(val), std::begin(local_val));


            static int n = N;
            Teuchos::reduceAll(*get(), Teuchos::REDUCE_SUM, n, &local_val[0], &val[0]);
            return ret_global;
        }

        template<typename T>
        inline T min(const T &val) const {

            T ret_global = 0.0;
            Teuchos::reduceAll(*get(), Teuchos::REDUCE_MIN, 1, &val, &ret_global);
            return ret_global;
        }

        template<typename T>
        inline T max(const T &val) const {

            T ret_global = 0.0;
            Teuchos::reduceAll(*get(), Teuchos::REDUCE_MAX, 1, &val, &ret_global);
            return ret_global;
        }

        inline static TrilinosCommunicator &get_default()
        {
            static TrilinosCommunicator instance_;
            return instance_;
        }

        TrilinosCommunicator split(const int color) const;

        TrilinosCommunicator(const CommPtr &comm) : comm_(comm) {}
        TrilinosCommunicator(const SelfCommunicator &comm);
        TrilinosCommunicator(const Communicator &comm);
        TrilinosCommunicator();

    private:
        CommPtr comm_;
    };
}


#endif //UTOPIA_TRILINOS_COMMUNICATOR_HPP
