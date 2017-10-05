#include "utopia_Log.hpp"
#include "utopia_Core.hpp"

#include <numeric>
#include <fstream>

#ifdef UTOPIA_LOG_ENABLED

namespace utopia {

    MeasurementId Measurement::generate_unique_id(){
        static long l = mpi_world_rank();
        return l += mpi_world_size();
    }

    Log::Log() {
        start_time_ = MPI_Wtime();
    }

    Log &Log::instance() {
        static Log instance;
        return instance;
    }

    void Log::save_collected_log() {
        std::map<std::string, std::vector<Measurement>> class_group;

        std::ofstream f_detail("log." + std::to_string(mpi_world_rank()) + ".csv");
        f_detail << std::fixed << std::setprecision(9) << "Time;Duration (s);Class" << std::endl;

        for (const auto &it : event_map_) {
            const Measurement &m = it.second;
            class_group[m.class_].push_back(m);

            double timestamp = m.start_time_ - start_time_;

            f_detail << timestamp << ';' << m.end_time_ - m.start_time_ << ';' << m.class_ << std::endl;
        }

        f_detail.close();

        std::ofstream f_summary("summary." + std::to_string(mpi_world_rank()) + ".csv");
        f_summary << std::fixed <<
            "Class;Total time spent (s);Count;Mean time (s);Std deviation (s);Relative std dev\n";

        for (const auto &it : class_group) {
            std::vector<double> values;
            values.reserve(it.second.size());
            std::transform(it.second.begin(), it.second.end(), std::back_inserter(values),
                [](const Measurement &m) {
                    return m.end_time_ - m.start_time_;
                }
            );

            double total_time = std::accumulate(values.begin(), values.end(), 0.0);
            double mean = total_time / values.size();

            double stddev = std::sqrt(std::accumulate(values.begin(), values.end(), 0.0,
                [mean](double v, double t) {
                    return v + (t - mean) * (t - mean);
                }) / values.size());

            f_summary << std::setprecision(9)
                << it.first << ';' << total_time << ';' << values.size()
                << ';' << mean << ';' << stddev << ';'
                << std::setprecision(2) << stddev / mean << std::endl;
        }

        f_summary.close();

		std::ofstream f_memory("memory." + std::to_string(mpi_world_rank()) + ".csv");
		f_memory << std::fixed << "Time (s);Memory used by pool\n";

		for (const auto &it : memory_log_) {
			double timestamp = it.first - start_time_;

			f_memory << timestamp << ';' << it.second << std::endl;
		}

		f_memory.close();
    }
}

#endif  //UTOPIA_LOG_ENABLED
