#include "utopia_TestRegistry.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Chrono.hpp"

#include <iostream>


namespace utopia {

    char TestRegistry::add_test_unit(const std::string &unit_name, TestRegistry::RunTest run_test)
    {
        instance().units_[unit_name] = run_test;
        return 0;
    }

    char TestRegistry::add_optional_test_unit(const std::string &unit_name, RunTest run_test)
    {
        instance().optional_units_[unit_name] = run_test;
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
        os << "select with: -test <sub-command>\n";
        os << "available sub-commands:\n";
        for(const auto &u : units_) {
            os << "\t" << u.first << "\n";
        }

        if(!optional_units_.empty()) {
            os << "\t------ optional -------" << std::endl;

            for(const auto &u : optional_units_) {
                os << "\t" << u.first << "\n";
            }
        }
        os << std::flush;
    }

    void TestRegistry::run_unit(const std::string &unit_name, RunTest run_test)
    {
        try {
            if(utopia::mpi_world_rank() == 0 && verbose()) {                                     \
                std::cout << "--------------------------------------------------------\n";       \
                std::cout << "begin:\t[" << unit_name << "] unit tests" << std::endl;  \
                std::cout << "--------------------------------------------------------\n";       \
            }     

            run_test();

            if(utopia::mpi_world_rank() == 0 && verbose()) {                                     \
                std::cout << "--------------------------------------------------------\n";       \
                std::cout << "end:\t[" << unit_name << "] unit tests" << std::endl;  \
                std::cout << "--------------------------------------------------------\n";       \
            }  

            } catch(std::exception &ex) {
                std::cerr << "[Failure] in " << unit_name << " " << ex.what() << std::endl;
                error_code_ = 1;
        } 
    }

    TestRegistry::~TestRegistry() {
        // if(mpi_world_rank() == 0) {
            std::cout << "[Status] test ran " << n_test_ran_ << std::endl;
        // }
    }

    int TestRegistry::run_all()
    {
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin testing]" << std::endl;
        }

        Chrono c;
        c.start();

        int n_runned = 0;
        for(std::map<std::string, RunTest>::const_iterator it = units_.begin(); it != units_.end(); ++it) {
            run_unit(it->first, it->second);
            ++n_runned;
        }

        if(mpi_world_rank() == 0) {
            std::cout << "[End testing" <<  ", Ran " << n_runned << " units]" << std::endl;
        }

        mpi_world_barrier();
        c.stop();
        if(mpi_world_rank() == 0) {
            std::cout << c << std::endl;
        }

        mpi_world_barrier();
        return error_code_;

    }

    int TestRegistry::run_aux(const std::map<std::string, RunTest> &units, const std::string &unit_name)
    {
        auto it = units.find(unit_name);
        if(it == units.end()) {
            return -1;
        }
        run_unit(it->first, it->second);
        return error_code_;
    }

    int TestRegistry::run(const std::string &unit_name)
    {
        int ret = run_aux(units_, unit_name);
        if(ret == -1) {
            ret = run_aux(optional_units_, unit_name);
        }

        if(ret == -1) {
            std::cerr << "[Error] no unit test with name " << unit_name << std::endl;
        }

        return ret;
    }

}

