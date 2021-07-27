#ifndef UTOPIA_TEST_ASSERT_HPP
#define UTOPIA_TEST_ASSERT_HPP

#include <sstream>

namespace utopia {
    class TestAssert {
    public:
        template <typename Left, typename Right>
        static void equal(const char *filename,
                          const int line,
                          const char *l_name,
                          const char *r_name,
                          const Left &l,
                          const Right &r) {
            if (l != r) {
                std::stringstream ss;
                ss << "Failure in expression: " << l_name << " == " << r_name << ".\n";
                ss << "l_name = " << l << '\n';
                ss << "r_name = " << r << '\n';
                failure(filename, line, ss.str());
            }
        }

        static void is_true(const char *filename, const int line, const char *expr, const bool &value);

        static void failure(const char *filename, const int line, const std::string &message);
    };

}  // namespace utopia

#define UTOPIA_TEST_EQ(_macro_left_, _macro_right_) \
    utopia::TestAssert::equal(__FILE__, __LINE__, #_macro_left_, #_macro_right_, _macro_left_, _macro_right_)

#define UTOPIA_TEST_TRUE(_macro_expr_) utopia::TestAssert::is_true(__FILE__, __LINE__, #_macro_expr_, _macro_expr_)

#endif  // UTOPIA_TEST_ASSERT_HPP
