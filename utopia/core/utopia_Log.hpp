#ifndef UTOPIA_UTOPIA_LOG_HPP
#define UTOPIA_UTOPIA_LOG_HPP

#include <iomanip>
#include <iostream>
#include <map>
#include <stack>
#include <chrono>

#include "utopia_Core.hpp"

#ifdef UTOPIA_LOG_ENABLED

namespace utopia {

    class Measurement;
    typedef long MeasurementId;

    class Log {
    public:
        template<class T>
        inline MeasurementId apply_begin(const Expression<T> &expr);

        inline void apply_end();

        static Log &instance();

        void save_collected_log();

    private:
        Log();
        std::chrono::high_resolution_clock::time_point start_time_;
        std::stack<MeasurementId> running_events_;
        std::map<MeasurementId, Measurement> event_map_;
    };


    class Measurement {
    public:
        template<class T>
        Measurement(const Expression<T> &expr) {
            id_ = generate_unique_id();
            size_ = size(expr.derived());
            class_ = expr.getClass();
        }

        inline MeasurementId get_id() const {
            return id_;
        }

        inline void begin() {
            start_time_ = std::chrono::high_resolution_clock::now();
        }

        inline void end() {
            end_time_ = std::chrono::high_resolution_clock::now();
        }

        friend void Log::save_collected_log();

    private:
        MeasurementId generate_unique_id();

        MeasurementId id_;
        std::string class_;
        Size size_;
        std::chrono::high_resolution_clock::time_point start_time_, end_time_;
    };

    template<class T>
    inline MeasurementId Log::apply_begin(const Expression<T> &expr) {
        Measurement m(expr);
        running_events_.push(m.get_id());
        event_map_.insert(std::make_pair(m.get_id(), m));
        event_map_.at(m.get_id()).begin();
        return m.get_id();
    }

    inline void Log::apply_end() {
        const MeasurementId &id = running_events_.top();
        event_map_.at(id).end();
        running_events_.pop();
    }

}

#define UTOPIA_LOG_BEGIN(expr)  utopia::Log::instance().apply_begin(expr)
#define UTOPIA_LOG_END(expr)    utopia::Log::instance().apply_end()


#else  //UTOPIA_LOG_ENABLED

#define UTOPIA_LOG_BEGIN(...)   0
#define UTOPIA_LOG_END(...)     0

#endif  //UTOPIA_LOG_ENABLED


#endif //UTOPIA_UTOPIA_LOG_HPP
