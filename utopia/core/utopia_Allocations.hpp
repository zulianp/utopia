#ifndef UTOPIA_ALLOCATIONS_HPP
#define UTOPIA_ALLOCATIONS_HPP

//FIXME removeme
#define ENABLE_NO_ALLOC_REGIONS

#ifdef ENABLE_NO_ALLOC_REGIONS
#define UTOPIA_NO_ALLOC_BEGIN(macro_name_) utopia::Allocations::instance().no_alloc_region_begin(macro_name_)
#define UTOPIA_NO_ALLOC_END() utopia::Allocations::instance().no_alloc_region_end()
#define UTOPIA_REPORT_ALLOC(macro_name_) utopia::Allocations::instance().report_alloc(macro_name_, __FILE__, __LINE__)

namespace utopia {

    class Allocations {
    public:
        using Counter = unsigned long long;

        static Allocations &instance()
        {
            static Allocations instance_;
            return instance_;
        }

        inline void report_alloc(const std::string &name, const std::string &file, int line_number)
        {
            ++count_;

            if(is_no_allocation_region_) {
                std::cerr << "[VIOLATION] allocation (" << name << ") in region (" << region_name_.top() << ") at " << file << ":" << line_number << std::endl;
                assert(false);
            }
        }

        inline void no_alloc_region_begin(const std::string &name)
        {
            ++is_no_allocation_region_;
            region_name_.push(name);
        }

        inline void no_alloc_region_end()
        {
            --is_no_allocation_region_;
            region_name_.pop();
        }

        ~Allocations()
        {
            std::cout << "[Status] total allocations " << count_ << std::endl;
        }

    private:
        int is_no_allocation_region_;
        Counter count_;
        std::stack<std::sting> region_name_;

        Allocations() : is_no_allocation_region_(0), count_(0)
        {}
    };

}

#else

#define UTOPIA_NO_ALLOC_BEGIN(...) ((void)0)
#define UTOPIA_NO_ALLOC_END() ((void)0)
#define UTOPIA_REPORT_ALLOC(...) ((void)0)

#endif

#endif //UTOPIA_ALLOCATIONS_HPP
