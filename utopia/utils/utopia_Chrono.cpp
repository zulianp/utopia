#include "utopia_Chrono.hpp"

#include <algorithm>
#include <iostream>

#ifdef WITH_MPI
#include "mpi.h"
#include "utopia_MPI.hpp"
#endif  // WITH_MPI

namespace utopia {
    void Chrono::start() {
        start_ = std::chrono::high_resolution_clock::now();
#ifdef WITH_MPI
        mpi_start_ = MPI_Wtime();
#endif  // WITH_MPI

#ifdef WIN32
        QueryPerformanceCounter(&realtime_start_);
#elif __APPLE__
        realtime_start_ = mach_absolute_time();
#else
        gettimeofday(&realtime_start_, NULL);
#endif  // WIN32
    }

#ifdef __APPLE__
    double Chrono::realtime_duration() const {
        uint64_t difference = realtime_end_ - realtime_start_;
        static double conversion = 0.0;

        if (conversion == 0.0) {
            mach_timebase_info_data_t info;
            kern_return_t err = mach_timebase_info(&info);

            if (err == 0) conversion = 1e-9 * (double)info.numer / (double)info.denom;
        }

        return conversion * (double)difference;
    }
#endif  //__APPLE__

    void Chrono::stop() {
        end_ = std::chrono::high_resolution_clock::now();
        duration_ = end_ - start_;

#ifdef WITH_MPI
        mpi_end_ = MPI_Wtime();
        mpi_duration_ = mpi_end_ - mpi_start_;
#endif  // WITH_MPI

#ifdef WIN32
        double start_time_ms = 0;
        double end_time_ms = 0;

        QueryPerformanceCounter(&realtime_end_);

        start_time_ms = realtime_start_.QuadPart * (1000000.0 / _frequency.QuadPart);
        end_time_ms = realtime_end_.QuadPart * (1000000.0 / _frequency.QuadPart);
        realtime_duration_ = end_time_ms - start_time_ms;
#elif __APPLE__
        realtime_end_ = mach_absolute_time();
        realtime_duration_ = realtime_duration();
#else
        double start_time_ms = 0;
        double end_time_ms = 0;

        gettimeofday(&realtime_end_, NULL);
        start_time_ms = (realtime_start_.tv_sec * 1000000.0) + realtime_start_.tv_usec;
        end_time_ms = (realtime_end_.tv_sec * 1000000.0) + realtime_end_.tv_usec;
        realtime_duration_ = (end_time_ms - start_time_ms) * 1e-6;
#endif  // WIN32
    }

    void Chrono::describe(std::ostream &os) const {
        os << "[Time elapsed] clock: " << std::to_string(duration_.count()) << " milliseconds,\t";
        os << " " << realtime_duration_ << " seconds.\t";

#ifdef WITH_MPI
        if (mpi_world_rank() == 0) {
            os << "MPI_Wtime: " << mpi_duration_ << " seconds.\n";
        }
#endif  // WITH_MPI
    }

    Chrono &Chrono::operator+=(const Chrono &other) {
        using std::max;
        using std::min;

#ifdef WITH_MPI
        mpi_start_ = min(mpi_start_, other.mpi_start_);
        mpi_end_ = max(mpi_end_, other.mpi_end_);
        mpi_duration_ += other.mpi_duration_;
#endif  // WITH_MPI

#ifdef __APPLE__
        realtime_start_ = min(realtime_start_, other.realtime_start_);
        realtime_end_ = max(realtime_end_, other.realtime_end_);
        realtime_duration_ += other.realtime_duration_;
#endif  //__APPLE__

#ifdef WIN32
        static bool not_impl_msg = true;
        if (not_impl_msg) {
            std::cerr << "[Warning] Chrono &Chrono::operator+=(const Chrono &other) not implemented for windows"
                      << std::endl;
            not_impl_msg = false;
        }
#endif

        // FIXME add Linux variant
        UTOPIA_UNUSED(other);

        return *this;
    }

    void Chrono::rescale_duration(const double factor) {
#ifdef WITH_MPI
        mpi_duration_ *= factor;
#endif  // WITH_MPI

#ifdef __APPLE__
        realtime_duration_ *= factor;
#endif  //__APPLE__

#ifdef WIN32
        static bool not_impl_msg = true;
        if (not_impl_msg) {
            std::cerr << "[Warning] void Chrono::rescale_duration(const double factor) not implemented for windows"
                      << std::endl;
            not_impl_msg = false;
        }
#endif

        // FIXME add Linux variant
        UTOPIA_UNUSED(factor);
    }

}  // namespace utopia
