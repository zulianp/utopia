#ifndef UTOPIA_CHRONO_HPP
#define UTOPIA_CHRONO_HPP

#include <chrono>
#include <ctime>
#include <ostream>

#ifdef WIN32   // Windows system specific
#include <windows.h>
#elif __APPLE__ // Unix based system specific
#include <mach/mach_time.h> // for mach_absolute_time
#else
#include <sys/time.h>
#endif

#include "utopia_Base.hpp"

namespace utopia {
	class Chrono {
	public:
		void start();
		void stop();
		void describe(std::ostream &os) const;
		
		inline double get_seconds() const
		{
			return realtime_duration_;
		}

		inline friend std::ostream & operator<<(std::ostream &os, const Chrono &c)
		{
			c.describe(os);
			return os;
		}

	private:
		typedef std::chrono::high_resolution_clock::time_point TimePoint;
		typedef std::chrono::duration<double, std::milli> DurationMillis;

		TimePoint start_, end_;
		DurationMillis duration_;

#ifdef WITH_MPI
		double mpi_start_, mpi_end_, mpi_duration_;
#endif //WITH_MPI		

#ifdef WIN32
		LARGE_INTEGER _frequency;      
		LARGE_INTEGER realtime_start_;     
		LARGE_INTEGER realtime_end_;       
#elif __APPLE__
		uint64_t realtime_start_;           
		uint64_t realtime_end_;  

		double realtime_duration() const;
#else
		timeval realtime_start_;           
		timeval realtime_end_;             
#endif

		double realtime_duration_;
	};
}

#endif //UTOPIA_CHRONO_HPP
