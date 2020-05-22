#include "utopia_Base.hpp"
#include "utopia_Instance.hpp"

#include <iostream>

#ifdef NDEBUG

namespace utopia {
    void test_check_assertion(const bool expr,
                              const std::string &filename,
                              const int line,
                              const std::string &expr_string) {
        if (!expr) {
            utopia::Utopia::instance().set_exit_code(1);
            std::cerr << filename << ": " << line << "\ntest failure: " << expr_string << std::endl;
#ifndef NDEBUG
            abort();
#endif
        }
    }
}  // namespace utopia

#endif  // NDEBUG
