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
        start_time_ = std::chrono::high_resolution_clock::now();
    }

    Log &Log::instance() {
        static Log instance;
        return instance;
    }

    void Log::save_collected_log() {
        std::map<std::string, std::vector<Measurement>> class_group;

        std::ofstream f_detail("log." + std::to_string(mpi_world_rank()) + ".csv");
        f_detail << "Time;Duration (us);Class" << std::endl;

        for (auto it = event_map_.cbegin(); it != event_map_.cend(); ++it) {
            const Measurement &m = it->second;
            class_group[m.class_].push_back(m);

            double timestamp = std::chrono::duration<double>(
                m.start_time_ - start_time_).count();

            f_detail << timestamp << ';' << std::chrono::high_resolution_clock::duration(
                m.end_time_ - m.start_time_).count() << ';' << m.class_ << std::endl;
        }

        f_detail.close();

        std::ofstream f_summary("summary." + std::to_string(mpi_world_rank()) + ".csv");
        f_summary << "Class;Total time spent (us);Count;Average time per operation (us)\n";

        for (auto it = class_group.cbegin(); it != class_group.cend(); ++it) {
            long total_time = std::accumulate(it->second.begin(), it->second.end(), 0,
                [](long v, const Measurement &m) {
                    return v +
                        std::chrono::high_resolution_clock::duration(
                        m.end_time_ - m.start_time_).count();
                });

            f_summary << it->first << ';' << total_time << ';' << it->second.size()
                << ';' << total_time / it->second.size() << std::endl;
        }

        f_summary.close();
    }
}

#endif  //UTOPIA_LOG_ENABLED
