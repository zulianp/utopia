#ifndef UTOPIA_BENCHMARK_HPP
#define UTOPIA_BENCHMARK_HPP

#include "utopia_Chrono.hpp"
#include <map>
#include <string>
#include <iostream>

namespace utopia {
	
	class Benchmark {
	public:
		Benchmark()
		: initialized_(false)
		{}

		virtual ~Benchmark() {}
		virtual void initialize() = 0;
		virtual std::string name() = 0;

		void benchmark_begin()
		{
			chrono_all_.start();
		}

		void benchmark_end()
		{
			chrono_all_.stop();

			if(mpi_world_rank() == 0) {
				std::cout << "---------------------\n";
				std::cout << name() << "\n";
				std::cout << "---------------------\n";
				to_csv(std::cout);
				std::cout << "---------------------\n";
			}
		}


		inline void clear()
		{
			current_ = "";
			measurements_.clear();
		}

		void register_experiment(const std::string &name, const std::function<void()> &fun)
		{
			experiments_[name] = fun;
		}

		void run()
		{
			if(!initialized_) { initialize(); }

			benchmark_begin();

			for(auto &e : experiments_) {
				begin_experiment(e.first);
				e.second();
				end_experiment();
			}

			benchmark_end();
		}

		inline void to_csv(std::ostream &os) const
		{
			os << "name,seconds\n";

			for(auto &m : measurements_)
			{
				os << m.first << "," << m.second << "\n";
			}

			os << "overall," << chrono_all_.get_seconds() << std::endl;
		}

	private:
		inline void begin_experiment(const std::string &name)
		{
			current_ = name;
			chrono_.start();
		}

		inline void end_experiment()
		{
			chrono_.stop();
			measurements_[current_] += chrono_.get_seconds();
		}


		Chrono chrono_all_;
		Chrono chrono_;
		std::string current_;
		std::map<std::string, double> measurements_;

		std::map<std::string, std::function<void()> > experiments_;
		bool initialized_;
	};


}

#endif //UTOPIA_BENCHMARK_HPP
