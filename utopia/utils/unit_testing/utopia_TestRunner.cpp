#include "utopia_TestRunner.hpp"
#include "utopia_Base.hpp"
#include "utopia_TestRegistry.hpp"

#include <iostream>

namespace utopia {

    TestRunner::TestRunner() = default;

    TestRunner::~TestRunner() = default;

    void TestRunner::verbose(const bool val) { TestRegistry::instance().verbose(val); }

    int TestRunner::run(int argc, char **argv) const {
        bool run_tests = true;

        std::vector<std::string> tests;
        for (int i = 1; i < argc; i++) {
            if (argv[i] == std::string("-test")) {
                if (++i >= argc) {
                    break;
                }
                tests.emplace_back(argv[i]);
            }

            else if (argv[i] == std::string("-list")) {
                this->describe();
                run_tests = false;
            }
        }

        if (run_tests) {
            if (tests.empty()) {
                return TestRegistry::instance().run_all();
            } else {
                return this->run(tests);
            }
        }

        return 0;
    }  // namespace utopia

    int TestRunner::run(const std::vector<std::string> &tests) const {
        auto &tr = TestRegistry::instance();
        int error_code = 0;
        for (const auto &t : tests) {
            int temp = 0;
            if ((temp = tr.run(t)) != 0) {
                error_code = temp;
            }
        }

        return error_code;
    }

    void TestRunner::describe(std::ostream &os) const { TestRegistry::instance().describe(os); }

}  // namespace utopia
