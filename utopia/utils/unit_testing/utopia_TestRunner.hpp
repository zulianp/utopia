#ifndef UTOPIA_TEST_RUNNER_HPP
#define UTOPIA_TEST_RUNNER_HPP

namespace utopia {

    class TestRunner {
    public:
        TestRunner();
        virtual ~TestRunner();
        int run(int argc, char **argv) const;
    };
}

#endif // UTOPIA_TEST_RUNNER_HPP
