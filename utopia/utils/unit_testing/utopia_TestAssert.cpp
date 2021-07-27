#include "utopia_TestAssert.hpp"

#include "utopia_IOStream.hpp"

namespace utopia {

    void TestAssert::is_true(const char *filename, const int line, const char *expr, const bool &value) {
        if (!value) {
            std::stringstream ss;
            ss << "Failure in expression: " << expr << ", expected true.\n";
            failure(filename, line, ss.str());
        }
    }

    void TestAssert::failure(const char *filename, const int line, const std::string &message) {
        utopia::err() << "At " << filename << ":" << line << ' ' << message << '\n';
    }

}  // namespace utopia
