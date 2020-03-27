#include "utopia_Tracer.hpp"
#include "utopia_Core.hpp"

#include <numeric>
#include <fstream>

#ifdef UTOPIA_TRACE_ENABLED

namespace utopia {

    MeasurementId Measurement::generate_unique_id(){
        static long l = mpi_world_rank();
        return l += mpi_world_size();
    }

    Tracer::Tracer() : full_trace_(false) {
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    Tracer &Tracer::instance() {
        static Tracer instance;
        return instance;
    }

    void Tracer::save_collected_log() {

        if(full_trace_) {
            std::map<std::string, std::vector<Measurement>> class_group;

            std::ofstream f_detail("log." + std::to_string(mpi_world_rank()) + ".csv");
            f_detail << std::fixed << std::setprecision(9) << "Time;Duration (s);Class" << std::endl;

            for (auto it = event_map_.cbegin(); it != event_map_.cend(); ++it) {
                const Measurement &m = it->second;
                class_group[m.class_].push_back(m);

                double timestamp = std::chrono::duration<double>(
                    m.start_time_ - start_time_).count();

                f_detail << timestamp << ';' << std::chrono::duration<double>(
                    m.end_time_ - m.start_time_).count() << ';' << m.class_ << std::endl;
            }

            f_detail.close();

            std::ofstream f_summary("summary." + std::to_string(mpi_world_rank()) + ".csv");
            std::ofstream f_problematic("problematic." + std::to_string(mpi_world_rank()) + ".csv");

            {
                std::stringstream ss;
                ss << std::fixed <<
                    "Class;Total time spent (s);Count;Allocations;Mean time (s);Std deviation (s);Relative std dev\n";

                std::string str = ss.str();
                f_summary << str;
                f_problematic << str;
            }



            for (auto it = class_group.cbegin(); it != class_group.cend(); ++it) {
                std::vector<double> values;

                values.reserve(it->second.size());

                std::transform(it->second.begin(), it->second.end(), std::back_inserter(values),
                    [](const Measurement &m) {
                        return std::chrono::duration<double>(m.end_time_ - m.start_time_).count();
                    }
                );

                auto n_allocs = std::accumulate(it->second.begin(), it->second.end(), Allocations::Counter(0),
                    [](const Allocations::Counter &val, const Measurement &m) -> Allocations::Counter {
                        return m.count_allocs_ + val;
                    }
                );

                double total_time = std::accumulate(values.begin(), values.end(), 0.0);
                double mean = total_time / values.size();

                double stddev = std::sqrt(std::accumulate(values.begin(), values.end(), 0.0,
                    [mean](double v, double t) {
                        return v + (t - mean) * (t - mean);
                    }) / values.size());

                std::stringstream ss;

                ss << std::setprecision(9)
                    << it->first << ';' << total_time << ';' << values.size() << ";" << n_allocs
                    << ';' << mean << ';' << stddev << ';'
                    << std::setprecision(2) << stddev / mean << std::endl;


                std::string str = ss.str();
                f_summary << str;

                if(n_allocs > values.size()) {
                    f_problematic << str;
                }

                ss.clear();
            }

            f_summary.close();
            f_problematic.close();
        } else {
            std::ofstream f_summary("summary." + std::to_string(mpi_world_rank()) + ".csv");

            const int rank = mpi_world_rank();
            if(rank == 0) {
                f_summary << "Class;MPI rank;";
                TraceSummary::header(f_summary);
            }

            for(auto it = summary_.cbegin(); it != summary_.cend(); ++it) {
                f_summary << it->first << ";" << rank << ";";
                it->second.describe(f_summary);
            }

            f_summary.close();
        }
    }
}

#endif  //UTOPIA_TRACE_ENABLED
