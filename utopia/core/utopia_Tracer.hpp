#ifndef UTOPIA_TRACER_HPP
#define UTOPIA_TRACER_HPP

#include <iomanip>
#include <iostream>
#include <map>
#include <stack>
#include <chrono>
#include <csignal>
#include <string>

#include "utopia_Base.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Allocations.hpp"

#ifdef UTOPIA_TRACE_ENABLED

namespace utopia {

    class Measurement;
    typedef long MeasurementId;

    class Interceptor {
    public:

        void interrupt_on_intercept(const bool val)
        {
            interrupt_on_intercept_ = val;
        }

        template<class T>
        void intercept(const Expression<T> &expr) 
        {
            if(interrupt_on_intercept_ && !expr_.empty() && expr.get_class() == expr_)  {
                // abort();
                std::raise(SIGINT);

            }
        }

        void intercept(const std::string &name) 
        {
            if(interrupt_on_intercept_ && !expr_.empty() && name == expr_)  {
                // abort();
                std::raise(SIGINT);

            }
        }

        inline void expr(const std::string &expr) {
            expr_ = expr;
        }

        Interceptor() : expr_(), interrupt_on_intercept_(false) {}

    private:

        std::string expr_;
        bool interrupt_on_intercept_;
        
    };

    class Tracer {
    public:
        template<class T>
        MeasurementId apply_begin(const Expression<T> &expr);
        MeasurementId region_begin(const std::string &region_name);//, const std::string &file, int line);

        template<class T>
        MeasurementId apply_begin_specialized(const Expression<T> &expr);
        inline void apply_end();
        inline void region_end();
        static Tracer &instance();

        void save_collected_log();
        inline Interceptor &interceptor() { return interceptor_; }

    private:
        Tracer();
        std::chrono::high_resolution_clock::time_point start_time_;
        std::stack<MeasurementId> running_events_;
        std::map<MeasurementId, Measurement> event_map_;

        Interceptor interceptor_;
    };


    class Measurement {
    public:
        // Calls to size(expr.derived()) can cause compilation errors even if we check if there is a
        // version of size that takes T, because that can be a recusive call that can fail later
        // (example: Wrapper<PetscSerialSparseMatrix>). Size logging has been removed.
        template<class T>
        Measurement(const Expression<T> &expr, const std::string &prefix = "") {
            id_ = generate_unique_id();
            class_ = prefix + expr.get_class();
        }

        inline Measurement(const std::string &region_name) {
            id_ = generate_unique_id();
            class_ = region_name;
        }

        inline MeasurementId get_id() const {
            return id_;
        }

        inline void begin() {
            start_time_ = std::chrono::high_resolution_clock::now();
            count_allocs_ = Allocations::instance().count();
        }

        inline void end() {
            end_time_ = std::chrono::high_resolution_clock::now();
            count_allocs_ = Allocations::instance().count() - count_allocs_;
        }

        friend void Tracer::save_collected_log();

    private:
        MeasurementId generate_unique_id();

        MeasurementId id_;
        std::string class_;
        std::chrono::high_resolution_clock::time_point start_time_, end_time_;
        Allocations::Counter count_allocs_;
    };

    template<class T>
    inline MeasurementId Tracer::apply_begin(const Expression<T> &expr) {
        interceptor().intercept(expr);

        Measurement m(expr);
        running_events_.push(m.get_id());
        event_map_.insert(std::make_pair(m.get_id(), m));
        event_map_.at(m.get_id()).begin();
        return m.get_id();
    }

    template<class T>
    inline MeasurementId Tracer::apply_begin_specialized(const Expression<T> &expr) {
        interceptor().intercept(expr);

        Measurement m(expr, "specialized_");
        running_events_.push(m.get_id());
        event_map_.insert(std::make_pair(m.get_id(), m));
        event_map_.at(m.get_id()).begin();
        return m.get_id();
    }

    inline MeasurementId Tracer::region_begin(const std::string &region_name)//, const std::string &file, int line)
    {
       interceptor().intercept(region_name); 

       Measurement m(region_name);// + "  [" + file + ":" + std::to_string(line) + "]");
       running_events_.push(m.get_id());
       event_map_.insert(std::make_pair(m.get_id(), m));
       event_map_.at(m.get_id()).begin();
       return m.get_id();
    }

    inline void Tracer::apply_end() {
        const MeasurementId &id = running_events_.top();
        event_map_.at(id).end();
        running_events_.pop();
    }

    inline void Tracer::region_end() {
       apply_end();
    }

}

#define UTOPIA_TRACE_BEGIN(macro_expr_)  utopia::Tracer::instance().apply_begin(macro_expr_)
#define UTOPIA_TRACE_END(macro_expr_)    utopia::Tracer::instance().apply_end()

#define UTOPIA_TRACE_REGION_BEGIN(macro_region_name_)  utopia::Tracer::instance().region_begin(macro_region_name_);//, __FILE__, __LINE__)
#define UTOPIA_TRACE_REGION_END(macro_region_name_)    utopia::Tracer::instance().region_end()

#define UTOPIA_TRACE_BEGIN_SPECIALIZED(macro_expr_)  utopia::Tracer::instance().apply_begin_specialized(macro_expr_)
#define UTOPIA_TRACE_END_SPECIALIZED(macro_expr_)    utopia::Tracer::instance().apply_end()


#else  //UTOPIA_TRACE_ENABLED

#define UTOPIA_TRACE_BEGIN(...)   ((void)0)
#define UTOPIA_TRACE_END(...)     ((void)0)

#define UTOPIA_TRACE_REGION_BEGIN(...) ((void)0)
#define UTOPIA_TRACE_REGION_END(...) ((void)0)

#define UTOPIA_TRACE_BEGIN_SPECIALIZED(...) ((void)0)
#define UTOPIA_TRACE_END_SPECIALIZED(...) ((void)0)

#endif  //UTOPIA_TRACE_ENABLED


#endif //UTOPIA_TRACER_HPP
