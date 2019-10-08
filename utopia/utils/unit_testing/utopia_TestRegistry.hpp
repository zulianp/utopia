#ifndef UTOPIA_TEST_REGISTRY_H
#define UTOPIA_TEST_REGISTRY_H 

#include <map>
#include <string>
#include <functional>
#include <iostream>

namespace utopia {

    class TestRegistry {
    public:
        using Count = long;

        typedef void (*RunTest)();

        static char add_test_unit(const std::string &unit_name, RunTest run_test);
        static char add_optional_test_unit(const std::string &unit_name, RunTest run_test);

        static TestRegistry &instance();
        int run(const std::string &unit_name);
        int run_all();
        void describe(std::ostream &os = std::cout) const;
        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) { verbose_ = val; }
        inline static void test_ran() { ++(instance().n_test_ran_); }

    private:
        TestRegistry() : verbose_(false), error_code_(0), n_test_ran_(0) {}
        ~TestRegistry();
        
        std::map<std::string, RunTest> units_;
        std::map<std::string, RunTest> optional_units_;
        bool verbose_;
        int error_code_ = 0;
        Count n_test_ran_;

        void run_unit(const std::string &unit_name, RunTest run_test);
        int run_aux(const std::map<std::string, RunTest> &units, const std::string &unit_name);
    };

}

#endif //UTOPIA_TEST_REGISTRY_H

