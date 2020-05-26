#include "utopia_Base.hpp"
#include "utopia_Instance.hpp"

#include <cmath>
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

    // template <class Left, class Right, typename Scalar>
    // void test_assert_approxeq(const Left &left,
    //                           const Right &right,
    //                           const Scalar &tol,
    //                           const std::string &filename,
    //                           const int line,
    //                           const std::string &expr_string) {
    //     if (std::abs(left - right) > tol) {
    //         std::cerr << "assertion failure: " << expr_string << std::endl;
    //         std::cerr << "at " << filename << ":" << line << std::endl;

    //         std::cerr << "====================================================\n";

    //         std::cerr << "Left: ";
    //         disp(left);

    //         std::cerr << "====================================================\n";

    //         std::cerr << "Right: ";
    //         disp(right);

    //         std::cerr << "====================================================\n";

    //         abort();
    //     }
    // }

}  // namespace utopia

#endif  // NDEBUG
