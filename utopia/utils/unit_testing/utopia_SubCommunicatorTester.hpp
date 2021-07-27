#ifndef UTOPIA_SUB_COMMUNICATOR_TESTER_HPP
#define UTOPIA_SUB_COMMUNICATOR_TESTER_HPP

#include "utopia_IOStream.hpp"
#include "utopia_UnitTest.hpp"

namespace utopia {

    template <class Comm>
    class SubCommunicatorTester {
    public:
        virtual ~SubCommunicatorTester() = default;

        void run(UnitTest<Comm> &test, const Comm &comm, const bool verbose = false) {
            const int comm_size = comm.size();
            const bool print_info = verbose && comm.rank() == 0;

            if (print_info) {
                utopia::out() << "--> [comm size = " << comm_size << "]" << std::endl;
            }

            // 1) we run on comm
            test.set_comm(comm);
            test.set_up();
            test.run();
            test.tear_down();

            // 1) we create sub-comms
            for (auto &p : parititions_sizes_) {
                if (p < comm_size) {
                    Comm sub_comm = comm.split(comm.rank() / p);

                    if (print_info) {
                        utopia::out() << "--> [sub comm(0) size: = " << sub_comm.size() << "]" << std::endl;
                    }

                    test.set_comm(sub_comm);
                    test.set_up();
                    test.run();
                    test.tear_down();
                } else {
                    break;
                }
            }
        }

        SubCommunicatorTester() : parititions_sizes_{{1, 2, 4}} {}

    private:
        std::vector<int> parititions_sizes_;
    };
}  // namespace utopia

#endif  // UTOPIA_SUB_COMMUNICATOR_TESTER_HPP
