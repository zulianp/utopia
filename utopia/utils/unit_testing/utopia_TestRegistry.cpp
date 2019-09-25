#include "utopia_TestRegistry.hpp"

#include <iostream>


namespace utopia {

    char TestRegistry::add_test_unit(const std::string &unit_name, TestRegistry::RunTest run_test)
    {
        instance().units_[unit_name] = run_test;
        return 0;
    }

    TestRegistry &TestRegistry::instance()
    {
        static TestRegistry instance_;
        return instance_;
    }

    void TestRegistry::describe(std::ostream &os) const
    {
        os << "Number of tests units: " << units_.size() << std::endl;
    }

    int TestRegistry::run_all()
    {
        int error_code = 0;
        for(std::map<std::string, RunTest>::const_iterator it = units_.begin(); it != units_.end(); ++it) {
            try {
                it->second();
            } catch(std::exception &ex) {
                std::cerr << "[Failure] in " << it->first << " " << ex.what() << std::endl;
                error_code = 1;
            }
        }

        return error_code;
    }

    static TestRegistry &registry = TestRegistry::instance();

}

