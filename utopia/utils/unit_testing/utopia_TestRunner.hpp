#ifndef UTOPIA_TEST_RUNNER_HPP
#define UTOPIA_TEST_RUNNER_HPP

#include <string>
#include <vector>
#include <iostream>

namespace utopia {

    class TestRunner {
    public:
        TestRunner();
        virtual ~TestRunner();
        int run(int argc, char **argv) const;
        int run(const std::vector<std::string> &tests) const;
        void describe(std::ostream &os = std::cout) const;
    };
}

#endif // UTOPIA_TEST_RUNNER_HPP
