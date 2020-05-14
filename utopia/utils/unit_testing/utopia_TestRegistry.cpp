#include "utopia_TestRegistry.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Chrono.hpp"

#include <iostream>

namespace utopia {

    TestRegistry::TestRegistry()
    {
        tests_.set_type("test");
        optional_tests_.set_type("test");
    }

    TestRegistry::~TestRegistry() = default;

    bool TestRegistry::verbose() const { return tests_.verbose(); }
    void TestRegistry::verbose(const bool val)
    {
        tests_.verbose(val);
        optional_tests_.verbose(val);
    }

    char TestRegistry::add_test_unit(const std::string &unit_name, TestRegistry::RunTest run_test)
    {
        instance().tests_.add_action(unit_name, run_test);
        return 0;
    }

    char TestRegistry::add_optional_test_unit(const std::string &unit_name, RunTest run_test)
    {
        instance().optional_tests_.add_action(unit_name, run_test);
        return 0;
    }

    TestRegistry &TestRegistry::instance()
    {
        static TestRegistry instance_;
        return instance_;
    }

    void TestRegistry::describe(std::ostream &os) const
    {
        os << "Number of tests units: " << tests_.size() << std::endl;
        os << "select with: -test <sub-command>\n";
        os << "available sub-commands:\n";
        tests_.describe(os);

        if(!optional_tests_.empty()) {
            os << "\t------ optional -------" << std::endl;
           optional_tests_.describe(os);
        }

        os << std::flush;
    }

    int TestRegistry::run_all()
    {
       return tests_.apply_all();

    }

    int TestRegistry::run(const std::string &unit_name)
    {
        int ret = tests_.apply(unit_name);
        if(ret == -1) {
            ret = optional_tests_.apply(unit_name);
        }

        if(ret == -1) {
            std::cerr << "[Error] no unit test with name " << unit_name << std::endl;
        }

        return ret;
    }

}

