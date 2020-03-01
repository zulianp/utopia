#ifndef UTOPIA_MPI_TIME_STATISTICS
#define UTOPIA_MPI_TIME_STATISTICS

#include <string>
#include <iostream>
#include <map>

#include "utopia_Chrono.hpp"
#include "utopia_Communicator.hpp"

namespace utopia {

    class MPITimeStatistics {
    public:

        void start()
        {
            comm_.barrier();
            c_.start();
        }

        void stop()
        {
            comm_.barrier();
            c_.stop();
        }

        void stop_and_collect(const std::string &collect_name)
        {
            stop();
            times_[std::to_string(counter_++) + ") " + collect_name] = c_.get_seconds();
        }

        void stop_collect_and_restart(const std::string &collect_name)
        {
            stop();
            times_[std::to_string(counter_++) + ") " + collect_name] = c_.get_seconds();
            start();
        }

        void add_stat(const std::string &name, const double &time)
        {
            times_[name] += time;
        }

        void describe(std::ostream &os) const
        {
            comm_.barrier();
            if(comm_.rank() == 0) {
                os << "------------------------------------------\n";
                os << "Timings (seconds), ";
                os << "comm_size : " << comm_.size() << "\n";
                for(auto p : times_) {
                    os << p.first << " : " << p.second << "\n";
                }

                os << "------------------------------------------\n";
            }
            comm_.barrier();
        }

        MPITimeStatistics(Communicator &comm) : comm_(comm), counter_(0) {}

    private:
        Communicator &comm_;
        Chrono c_;
        std::map<std::string, double> times_;
        int counter_;
    };

}

#endif //UTOPIA_MPI_TIME_STATISTICS
