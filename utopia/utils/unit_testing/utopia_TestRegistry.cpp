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

    TestRegistry &TestRegistry::instance()
    {
        static TestRegistry instance_;
        return instance_;
    }

    void TestRegistry::describe(std::ostream &os) const
    {
        os << "Number of tests units: " << units_.size() << std::endl;
        os << "select with: -test <sub-command>";
        os << "sub-commands:\n";
        for(const auto &u : units_) {
            os << "\t" << u.first << "\n";
        }
        os << std::flush;
    }

    int TestRegistry::run_all()
    {
        if(mpi_world_rank() == 0) {
            std::cout << "[Begin testing]" << std::endl;
        }

        Chrono c;
        c.start();

        int error_code = 0;
        for(std::map<std::string, RunTest>::const_iterator it = units_.begin(); it != units_.end(); ++it) {
            try {

                if(utopia::mpi_world_rank() == 0 && verbose()) {                                     \
                    std::cout << "--------------------------------------------------------\n";       \
                    std::cout << "begin:\t[" << (it->first) << "] unit tests" << std::endl;  \
                    std::cout << "--------------------------------------------------------\n";       \
                }     

                it->second();

                if(utopia::mpi_world_rank() == 0 && verbose()) {                                     \
                    std::cout << "--------------------------------------------------------\n";       \
                    std::cout << "end:\t[" << (it->first) << "] unit tests" << std::endl;  \
                    std::cout << "--------------------------------------------------------\n";       \
                }    

            } catch(std::exception &ex) {
                std::cerr << "[Failure] in " << it->first << " " << ex.what() << std::endl;
                error_code = 1;
            }
        }

        if(mpi_world_rank() == 0) {
            std::cout << "[End testing]" << std::endl;
        }

        mpi_world_barrier();
        c.stop();
        if(mpi_world_rank() == 0) {
            std::cout << c << std::endl;
        }

        mpi_world_barrier();

        return error_code;
    }

    int TestRegistry::run(const std::string &unit_name)
    {
        auto it = units_.find(unit_name);
        if(it == units_.end()) {
            std::cerr << "[Error] could not find test with name " << unit_name << std::endl;
            return -1;
        }

        try {
            it->second();
        } catch(std::exception &ex) {
            std::cerr << "[Failure] in " << it->first << " " << ex.what() << std::endl;
            return 1;
        }

        return 0;
    }

}

