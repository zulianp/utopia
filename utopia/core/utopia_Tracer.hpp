#ifndef UTOPIA_TRACER_HPP
#define UTOPIA_TRACER_HPP

#include <chrono>
#include <csignal>
#include <iomanip>
#include <iostream>
#include <map>
#include <stack>
#include <string>

#include "utopia_Allocations.hpp"
#include "utopia_Base.hpp"
#include "utopia_Expression.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_MPI.hpp"
#include "utopia_Utils.hpp"

#ifdef UTOPIA_TRACE_ENABLED

namespace utopia {

    class Measurement;
    class TraceSummary;

    typedef long MeasurementId;

    class Interceptor {
    public:
        void interrupt_on_intercept(const bool val) { interrupt_on_intercept_ = val; }

        template <class T>
        void intercept(const Expression<T> &expr) {
            if (interrupt_on_intercept_ && !expr_.empty() && expr.get_class() == expr_) {
                // abort();
                std::raise(SIGINT);
            }
        }

        void intercept(const std::string &name) {
            if (interrupt_on_intercept_ && !expr_.empty() && name == expr_) {
                // abort();
                std::raise(SIGINT);
            }
        }

        inline void expr(const std::string &expr) { expr_ = expr; }

        Interceptor() : expr_(), interrupt_on_intercept_(false) {}

    private:
        std::string expr_;
        bool interrupt_on_intercept_;
    };

    class Tracer {
    public:
        template <class T>
        MeasurementId apply_begin(const Expression<T> &expr);
        MeasurementId region_begin(const std::string &region_name);  //, const std::string &file, int line);

        template <class T>
        MeasurementId apply_begin_specialized(const Expression<T> &expr);
        inline void apply_end();
        inline void region_end();
        static Tracer &instance();

        void save_collected_log();
        inline Interceptor &interceptor() { return interceptor_; }

        inline void full_trace(const bool val) { full_trace_ = val; }
        inline void only_root_export(const bool val) { only_root_export_ = val; }
        inline void log_regions(const bool val) { log_regions_ = val; }

    private:
        Tracer();
        std::chrono::high_resolution_clock::time_point start_time_;
        std::stack<MeasurementId> running_events_;
        std::map<MeasurementId, Measurement> event_map_;
        std::map<std::string, TraceSummary> summary_;

        Interceptor interceptor_;
        bool full_trace_;
        bool only_root_export_;
        bool log_regions_{false};
        int print_indent_level{0};
    };

    class Measurement {
    public:
        // Calls to size(expr.derived()) can cause compilation errors even if we check if there is a
        // version of size that takes T, because that can be a recusive call that can fail later
        // (example: Wrapper<PetscSerialSparseMatrix>). Size logging has been removed.
        template <class T>
        Measurement(const Expression<T> &expr, const std::string &prefix = "") {
            id_ = generate_unique_id();
            class_ = prefix + expr.get_class();
        }

        inline Measurement(const std::string &region_name) {
            id_ = generate_unique_id();
            class_ = region_name;
        }

        inline MeasurementId get_id() const { return id_; }

        inline const std::string &get_class() const { return class_; }

        inline void begin() {
            start_time_ = std::chrono::high_resolution_clock::now();
            count_allocs_ = Allocations::instance().count();
        }

        inline void end() {
            end_time_ = std::chrono::high_resolution_clock::now();
            count_allocs_ = Allocations::instance().count() - count_allocs_;
        }

        inline double seconds() const { return std::chrono::duration<double>(end_time_ - start_time_).count(); }

        friend void Tracer::save_collected_log();
        friend class TraceSummary;

    private:
        MeasurementId generate_unique_id();

        MeasurementId id_;
        std::string class_;
        std::chrono::high_resolution_clock::time_point start_time_, end_time_;
        Allocations::Counter count_allocs_;
    };

    class TraceSummary {
    public:
        inline TraceSummary &operator+=(const Measurement &m) {
            seconds_ += std::chrono::duration<double>(m.end_time_ - m.start_time_).count();
            count_allocs_ += m.count_allocs_;
            ++invokations_;

            return *this;
        }

        inline static void header(std::ostream &os) { os << "Total time spent (s);Count;Allocations;Mean time (s)\n"; }

        inline void describe(std::ostream &os) const {
            os << seconds_ << ";" << invokations_ << ";" << count_allocs_ << ";" << (seconds_ / invokations_) << "\n";
        }

        inline TraceSummary() : count_allocs_(0), invokations_(0), seconds_(0) {}

        inline double seconds() const { return seconds_; }

    private:
        Allocations::Counter count_allocs_;
        Allocations::Counter invokations_;
        double seconds_;
    };

    template <class T>
    inline MeasurementId Tracer::apply_begin(const Expression<T> &expr) {
        interceptor().intercept(expr);

        Measurement m(expr);

        if (log_regions_) {
            std::string indent(print_indent_level * 2, ' ');
            utopia::out() << "BEGIN[" << mpi_world_rank() << "]: " << m.get_class() << "\n";
            print_indent_level++;
        }

        running_events_.push(m.get_id());
        auto it = event_map_.insert(std::make_pair(m.get_id(), m));
        it.first->second.begin();
        return m.get_id();
    }

    template <class T>
    inline MeasurementId Tracer::apply_begin_specialized(const Expression<T> &expr) {
        interceptor().intercept(expr);

        Measurement m(expr, "specialized_");

        if (log_regions_) {
            std::string indent(print_indent_level * 2, ' ');
            utopia::out() << indent << "BEGIN[" << mpi_world_rank() << "]: " << m.get_class() << "\n";
            print_indent_level++;
        }

        running_events_.push(m.get_id());
        auto it = event_map_.insert(std::make_pair(m.get_id(), m));
        it.first->second.begin();
        return m.get_id();
    }

    inline MeasurementId Tracer::region_begin(const std::string &region_name)  //, const std::string &file, int line)
    {
        interceptor().intercept(region_name);

        Measurement m(region_name);  // + "  [" + file + ":" + std::to_string(line) + "]");

        if (log_regions_) {
            std::string indent(print_indent_level * 2, ' ');
            utopia::out() << indent << "BEGIN[" << mpi_world_rank() << "]: " << m.get_class() << "\n";
            print_indent_level++;
        }

        running_events_.push(m.get_id());
        auto it = event_map_.insert(std::make_pair(m.get_id(), m));
        it.first->second.begin();
        return m.get_id();
    }

    inline void Tracer::apply_end() {
        const MeasurementId &id = running_events_.top();
        auto it = event_map_.find(id);

        it->second.end();
        running_events_.pop();

        if (log_regions_) {
            print_indent_level--;
            std::string indent(print_indent_level * 2, ' ');
            utopia::out() << indent << "END[" << mpi_world_rank() << "]: " << it->second.get_class() << " ("
                          << it->second.seconds() << " seconds)\n";
        }

        if (!full_trace_) {
            summary_[it->second.get_class()] += it->second;
            event_map_.erase(it);
        }
    }

    inline void Tracer::region_end() { apply_end(); }

}  // namespace utopia

#ifdef UTOPIA_TRACE_EXPR_ENABLED
#define UTOPIA_TRACE_BEGIN(macro_expr_) utopia::Tracer::instance().apply_begin(macro_expr_)
#define UTOPIA_TRACE_END(macro_expr_) utopia::Tracer::instance().apply_end()
#define UTOPIA_TRACE_BEGIN_SPECIALIZED(macro_expr_) utopia::Tracer::instance().apply_begin_specialized(macro_expr_)
#define UTOPIA_TRACE_END_SPECIALIZED(macro_expr_) utopia::Tracer::instance().apply_end()
#else
#define UTOPIA_TRACE_BEGIN(...) ((void)0)
#define UTOPIA_TRACE_END(...) ((void)0)
#define UTOPIA_TRACE_BEGIN_SPECIALIZED(...) ((void)0)
#define UTOPIA_TRACE_END_SPECIALIZED(...) ((void)0)
#endif

#define UTOPIA_TRACE_SCOPE(macro_expr_)                   \
    utopia::Tracer::instance().region_begin(macro_expr_); \
    utopia::ExitScopeCallback utopia_trace_scope__([]() { utopia::Tracer::instance().region_end(); });

#define UTOPIA_TRACE_REGION_BEGIN(macro_region_name_) \
    utopia::Tracer::instance().region_begin(macro_region_name_);  //, __FILE__, __LINE__)
#define UTOPIA_TRACE_REGION_END(macro_region_name_) utopia::Tracer::instance().region_end()
#else  // UTOPIA_TRACE_ENABLED

#define UTOPIA_TRACE_SCOPE(macro_expr_) ((void)0)
#define UTOPIA_TRACE_BEGIN(...) ((void)0)
#define UTOPIA_TRACE_END(...) ((void)0)

#define UTOPIA_TRACE_REGION_BEGIN(...) ((void)0)
#define UTOPIA_TRACE_REGION_END(...) ((void)0)

#define UTOPIA_TRACE_BEGIN_SPECIALIZED(...) ((void)0)
#define UTOPIA_TRACE_END_SPECIALIZED(...) ((void)0)

#endif  // UTOPIA_TRACE_ENABLED

#endif  // UTOPIA_TRACER_HPP
