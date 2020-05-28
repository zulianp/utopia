#ifndef UTOPIA_PETSC_ASSERT_HPP
#define UTOPIA_PETSC_ASSERT_HPP

#include <petscsys.h>
#include "utopia_Assert.hpp"
#include "utopia_Base.hpp"

namespace utopia {

    template <class TensorLeft,
              class TensorRight,
              int Order = Traits<TensorLeft>::Order,
              int Backend = Traits<TensorLeft>::Backend>
    class AssertApproxEq<TensorLeft, TensorRight, Order, PETSC> {
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
                CHKERRQ(1);
            }

            return ok;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_ASSERT_HPP
