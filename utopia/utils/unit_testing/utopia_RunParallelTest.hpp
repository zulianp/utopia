#ifndef UTOPIA_RUN_PARALLEL_TEST_HPP
#define UTOPIA_RUN_PARALLEL_TEST_HPP

#include "utopia_SubCommunicatorTester.hpp"

namespace utopia {
    template <class Comm>
    void run_parallel_test(UnitTest<Comm> &test, bool verbose, const Comm &comm = Comm()) {
        test.print_backend_info();
        SubCommunicatorTester<Comm>().run(test, comm, verbose);
    }

    template <class Test>
    void run_parallel_test(bool verbose = false) {
        Test test;
        run_parallel_test(test, verbose);
    }

    template <class Test>
    void run_serial_test() {
        Test test;
        test.print_backend_info();
        test.run();
    }
}  // namespace utopia

#endif  // UTOPIA_RUN_PARALLEL_TEST_HPP