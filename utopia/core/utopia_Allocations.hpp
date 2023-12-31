#ifndef UTOPIA_ALLOCATIONS_HPP
#define UTOPIA_ALLOCATIONS_HPP

#include <iostream>
#include <stack>
#include <string>

#include "utopia_Base.hpp"
#include "utopia_IOStream.hpp"

#ifdef ENABLE_NO_ALLOC_REGIONS
#define UTOPIA_NO_ALLOC_BEGIN(macro_name_) utopia::Allocations::instance().no_alloc_region_begin(macro_name_)
#define UTOPIA_NO_ALLOC_END() utopia::Allocations::instance().no_alloc_region_end()
#define UTOPIA_REPORT_ALLOC(macro_name_) utopia::Allocations::instance().report_alloc(macro_name_, __FILE__, __LINE__)
#else
#define UTOPIA_NO_ALLOC_BEGIN(...) ((void)0)
#define UTOPIA_NO_ALLOC_END() ((void)0)
#define UTOPIA_REPORT_ALLOC(...) ((void)0)
#endif

namespace utopia {

    class Allocations final {
    public:
        using Counter = unsigned long long;

        inline static Allocations &instance() {
            static Allocations instance_;
            return instance_;
        }

        inline void report_alloc(const std::string &name, const std::string &file, int line_number) {
            ++count_;

            if (is_no_allocation_region_) {
                handle_violation(name, file, line_number);
            }
        }

        void handle_violation(const std::string &name, const std::string &file, int line_number) {
            ++n_violations_;

            if (abort_on_violation_ || verbose_) {
                std::cerr << "[VIOLATION] allocation (" << name << ") in region (" << region_name_.top() << ") at "
                          << file << ":" << line_number << std::endl;
            }

            if (abort_on_violation_) {
                assert(false);
                abort();
            }
        }

        inline void no_alloc_region_begin(const std::string &name) {
            ++is_no_allocation_region_;
            region_name_.push(name);
        }

        inline void no_alloc_region_end() {
            --is_no_allocation_region_;
            region_name_.pop();
        }

        inline ~Allocations() {
            if (!region_name_.empty()) {
                std::cerr << "[Error] incorrect regions are present in the code. number of open regions is "
                          << region_name_.size() << std::endl;
            }

            if (count_ > 0) {
                utopia::out() << "[Status] total allocations " << count_ << ", " << n_violations_ << " violations "
                              << std::endl;
            }
        }

        inline void abort_on_violation(const bool val) { abort_on_violation_ = val; }

        inline Counter count() const { return count_; }

        inline void verbose(const bool val) { verbose_ = val; }

    private:
        int is_no_allocation_region_{0};
        Counter count_{0};
        Counter n_violations_{0};
        std::stack<std::string> region_name_;
        bool abort_on_violation_{false};
        bool verbose_{true};

        inline Allocations() = default;
    };

}  // namespace utopia

#endif  // UTOPIA_ALLOCATIONS_HPP
