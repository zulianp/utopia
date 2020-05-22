#ifndef UTOPIA_ASSERT_HPP
#define UTOPIA_ASSERT_HPP

#include "utopia_Base.hpp"
#include "utopia_Operations.hpp"
#include "utopia_Traits.hpp"

#include <cassert>
#include <sstream>

namespace utopia {

    template <class TensorLeft,
              class TensorRight,
              int Order = Traits<TensorLeft>::Order,
              int Backend = Traits<TensorLeft>::Backend>
    class AssertApproxEq {
    public:
        template <typename Scalar>
        static bool apply(const TensorLeft &l,
                          const TensorRight &r,
                          const Scalar &tol,
                          const std::string &filename,
                          const int line) {
            bool ok = approxeq(l, r, tol);

            if (!ok) {
                auto &comm = l.comm();
                std::stringstream ss;

                ss << "assertion failure: "
                   << "left != right" << std::endl;
                ss << "at " << filename << ":" << line << std::endl;

                ss << "====================================================\n";

                ss << "Left: ";

                comm.synched_print(ss.str(), std::cerr);

                disp(l);

                ss.clear();

                ss << "====================================================\n";

                ss << "Right: ";

                comm.synched_print(ss.str(), std::cerr);

                disp(r);

                ss.clear();

                ss << "====================================================\n";
                comm.synched_print(ss.str(), std::cerr);
                // FIXME
                abort();
            }

            return ok;
        }
    };

    template <class Left, class Right, int Backend>
    class AssertApproxEq<Left, Right, 0, Backend> {
    public:
        using Scalar = decltype(Left() + Right());

        static bool apply(const Left &l,
                          const Right &r,
                          const Scalar &tol,
                          const std::string &filename,
                          const int line) {
            bool ok = approxeq(l, r, tol);

            if (!ok) {
                auto &out = std::cerr;

                out << "assertion failure: "
                    << "left != right" << std::endl;
                out << "at " << filename << ":" << line << std::endl;

                out << "\n====================================================\n";

                out << "Left:\n";

                out << l;

                out << "\n====================================================\n";

                out << "Right:\n";

                out << r;

                out << "\n====================================================\n";
                out << std::flush;
                // FIXME
                abort();
            }

            return ok;
        }
    };

    template <class TensorLeft, class TensorRight, typename Scalar>
    bool assert_approx_eq(const TensorLeft &l,
                          const TensorRight &r,
                          const Scalar &tol,
                          const std::string &filename,
                          const int line) {
        return AssertApproxEq<TensorLeft, TensorRight>::apply(l, r, tol, filename, line);
    }

}  // namespace utopia

#define utopia_test_asserteq(macro_left_, macro_right_, macro_tol_) \
    assert_approx_eq(macro_left_, macro_right_, macro_tol_, __FILE__, __LINE__)

#endif  // UTOPIA_ASSERT_HPP
