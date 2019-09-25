#ifndef UTOPIA_TEST_REGISTRY_H
#define UTOPIA_TEST_REGISTRY_H 

#include <map>
#include <string>
#include <functional>
#include <iostream>

namespace utopia {

    class TestRegistry {
    public:
        typedef void (*RunTest)();

        static char add_test_unit(const std::string &unit_name, RunTest run_test);

        static TestRegistry &instance();

        int run_all();
        void describe(std::ostream &os = std::cout) const;

    private:
        TestRegistry() {}
        std::map<std::string, RunTest> units_;
    };

}

#endif //UTOPIA_TEST_REGISTRY_H

