#ifndef UTOPIA_TEST_REGISTRY_H
#define UTOPIA_TEST_REGISTRY_H

#include <map>
#include <string>
#include <functional>
#include <iostream>

#include "utopia_ActionRegistry.hpp"

namespace utopia {

    class TestRegistry {
    public:
        using Count = long;

        using RunTest = void (*)();

        static char add_test_unit(const std::string &unit_name, RunTest run_test);
        static char add_optional_test_unit(const std::string &unit_name, RunTest run_test);

        static TestRegistry &instance();
        int run(const std::string &unit_name);
        int run_all();
        void describe(std::ostream &os = std::cout) const;
        bool verbose() const;
        void verbose(const bool val);
        static void test_ran();

    private:
        TestRegistry();
        ~TestRegistry();
        ActionRegistry tests_, optional_tests_;
    };

}

#endif //UTOPIA_TEST_REGISTRY_H

