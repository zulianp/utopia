#ifndef UTOPIA_COMMUNICATOR_HPP
#define UTOPIA_COMMUNICATOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Clonable.hpp"
#include "utopia_Traits.hpp"

#include "utopia_MPI_Operations.hpp"

#ifdef UTOPIA_WITH_MPI
#include <mpi.h>
#endif

namespace utopia {

    class Communicator : public Clonable {
    public:
        ~Communicator() override = default;
        virtual int rank() const = 0;
        virtual int size() const = 0;
        Communicator *clone() const override = 0;
        virtual void barrier() const = 0;

        virtual bool same(const Communicator &other) const { return size() == other.size(); }

#ifdef UTOPIA_WITH_MPI
        virtual MPI_Comm raw_comm() const = 0;
#endif

        virtual bool conjunction(const bool &val) const = 0;
        virtual bool disjunction(const bool &val) const = 0;

        inline bool is_root() const { return rank() == 0; }

        template <typename T>
        void root_print(const T &obj, std::ostream &os = std::cout) const {
            barrier();
            if (rank() == 0) {
                os << obj << std::endl;
            }
            barrier();
        }

        template <typename T>
        void synched_print(const T &obj, std::ostream &os = std::cout) const {
            const int n = size();
            const int r = rank();

            for (int i = 0; i < n; ++i) {
                barrier();

                if (i == r) {
                    os << "[" << r << "] ---------------------\n";
                    os << obj << std::endl;
                    os << std::flush;
                }
            }

            barrier();
        }

        // #ifdef UTOPIA_WITH_MPI
        // 		virtual MPI_Comm get() const = 0;
        // #endif //UTOPIA_WITH_MPI
    };

    /**
     * @brief usefull for non-mpi backends as a default object to return/provide
     */
    class SelfCommunicator final : public Communicator {
    public:
        int rank() const noexcept override { return 0; }
        int size() const noexcept override { return 1; }

        inline static SelfCommunicator world() { return SelfCommunicator(); }

        inline static SelfCommunicator self() { return SelfCommunicator(); }

        void barrier() const override {}

        SelfCommunicator *clone() const noexcept override { return new SelfCommunicator(); }

        inline static SelfCommunicator &get_default() {
            static SelfCommunicator instance_;
            return instance_;
        }

        inline bool conjunction(const bool &val) const override { return val; }

        inline bool disjunction(const bool &val) const override { return val; }

        template <typename T>
        inline T sum(const T &val) const noexcept {
            return val;
        }

        template <typename T>
        inline T min(const T &val) const noexcept {
            return val;
        }

        template <typename T>
        inline T max(const T &val) const noexcept {
            return val;
        }

#ifdef UTOPIA_WITH_MPI
        inline MPI_Comm get() const noexcept { return MPI_COMM_SELF; }

        inline MPI_Comm raw_comm() const override { return MPI_COMM_SELF; }

#endif  // UTOPIA_WITH_MPI
    };

    template <typename T>
    class MPIType {};

#ifdef UTOPIA_WITH_MPI
    template <>
    class MPIType<double> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_DOUBLE; }
    };

    template <>
    class MPIType<float> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_FLOAT; }
    };

    template <>
    class MPIType<long> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_LONG; }
    };

    template <>
    class MPIType<long long> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_LONG_LONG; }
    };

    template <>
    class MPIType<int> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_INT; }
    };

    template <>
    class MPIType<char> {
    public:
        inline static MPI_Datatype value() noexcept { return MPI_CHAR; }
    };

    class MPICommunicator : public Communicator {
    public:
        ~MPICommunicator() override = default;

        virtual MPI_Comm get() const = 0;

        /////////////////////////////////////////////

        bool conjunction(const bool &val) const override;
        bool disjunction(const bool &val) const override;

        inline void barrier() const override { MPI_Barrier(get()); }

        inline int rank() const override {
            int ret;
            MPI_Comm_rank(get(), &ret);
            return ret;
        }

        inline int size() const override {
            int ret;
            MPI_Comm_size(get(), &ret);
            return ret;
        }

        template <typename T>
        inline T sum(const T &val) const {
            T ret = val;
            MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_SUM, get());
            return ret;
        }

        template <typename T>
        inline void sum(const int n_values, T *inout) const {
            MPI_Allreduce(MPI_IN_PLACE, inout, n_values, MPIType<T>::value(), MPI_SUM, get());
        }

        template <typename T>
        inline void min(const int n_values, T *inout) const {
            MPI_Allreduce(MPI_IN_PLACE, inout, n_values, MPIType<T>::value(), MPI_MIN, get());
        }

        template <typename T>
        inline void max(const int n_values, T *inout) const {
            MPI_Allreduce(MPI_IN_PLACE, inout, n_values, MPIType<T>::value(), MPI_MAX, get());
        }

        template <typename T>
        inline T min(const T &val) const {
            T ret = val;
            MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_MIN, get());
            return ret;
        }

        template <typename T>
        inline T max(const T &val) const {
            T ret = val;
            MPI_Allreduce(MPI_IN_PLACE, &ret, 1, MPIType<T>::value(), MPI_MAX, get());
            return ret;
        }

        template <typename T>
        inline void exscan(const T *input, T *output, const int n, MPI_Op op) const {
            MPI_Exscan(input, output, n, MPIType<T>::value(), op, get());
        }

        template <typename T>
        inline void exscan_sum(const T *input, T *output, const int n) const {
            MPI_Exscan(input, output, n, MPIType<T>::value(), MPI_SUM, get());
        }

        template <class Op, typename T>
        inline T reduce(const Op & /*op*/, const T &val) const {
            T ret = val;
            MPI_Allreduce(MPI_IN_PLACE, &ret, 1, (MPIType<T>::value()), (MPIReduceOp<Op, HOMEMADE>::op()), get());
            return ret;
        }
    };

#endif  // UTOPIA_WITH_MPI

}  // namespace utopia

#endif  // UTOPIA_COMMUNICATOR_HPP
