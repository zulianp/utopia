#ifndef UTOPIA_BENCHMARK_HPP
#define UTOPIA_BENCHMARK_HPP

#include <iostream>
#include <map>
#include <string>
#include "utopia_Chrono.hpp"

namespace utopia {

    class Benchmark {
    public:
        Benchmark() = default;

        virtual ~Benchmark() = default;
        virtual void initialize() = 0;
        virtual std::string name() = 0;

        void benchmark_begin() { chrono_all_.start(); }

        void benchmark_end() {
            chrono_all_.stop();

            if (mpi_world_rank() == 0) {
                std::cout << "---------------------\n";
                std::cout << "[" << name() << "]"
                          << "\n";
                std::cout << "---------------------\n";
                to_csv(std::cout);
                std::cout << "---------------------\n";
            }
        }

        inline void clear() {
            current_ = "";
            measurements_.clear();
        }

        void register_experiment(const std::string &name, const std::function<void()> &fun) {
            experiments_[name] = fun;
        }

        void run() {
            if (!initialized_) {
                initialize();
            }

            benchmark_begin();

            for (auto &e : experiments_) {
                begin_experiment(e.first);
                e.second();
                end_experiment();
            }

            benchmark_end();
        }

        inline void to_csv(std::ostream &os) const {
            if (verbosity_level_ > 1) {
                os << "name,seconds\n";

                for (auto &m : measurements_) {
                    os << m.first << "," << m.second << "\n";
                }
            }

            os << "overall," << chrono_all_.get_seconds() << std::endl;
        }

        inline void set_verbosity_level(const int verbosity_level) { verbosity_level_ = verbosity_level; }

    private:
        inline void begin_experiment(const std::string &name) {
            current_ = name;
            chrono_.start();
        }

        inline void end_experiment() {
            chrono_.stop();
            measurements_[current_] += chrono_.get_seconds();
        }

        Chrono chrono_all_;
        Chrono chrono_;
        std::string current_;
        std::map<std::string, double> measurements_;

        std::map<std::string, std::function<void()> > experiments_;
        bool initialized_{false};
        int verbosity_level_{1};
    };

}  // namespace utopia

#endif  // UTOPIA_BENCHMARK_HPP
